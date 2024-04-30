import os
import os.path as osp
import argparse

from PIL import Image
import torch
from tqdm.auto import tqdm

import numpy as np
import matplotlib.pyplot as plt

from point_e.diffusion.configs import DIFFUSION_CONFIGS, diffusion_from_config
from point_e.diffusion.sampler import PointCloudSampler
from point_e.models.download import load_checkpoint
from point_e.models.configs import MODEL_CONFIGS, model_from_config
from point_e.util.plotting import plot_point_cloud

def main(args):
    dalle_folder = 'images/postprocessed_dalle' if args.cleaned else 'images/dalle'
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')

    print('creating base model...')
    base_name = args.pointe_name # use base300M or base1B for better results
    base_model = model_from_config(MODEL_CONFIGS[base_name], device)
    base_model.eval()
    base_diffusion = diffusion_from_config(DIFFUSION_CONFIGS[base_name])

    print('creating upsample model...')
    upsampler_model = model_from_config(MODEL_CONFIGS['upsample'], device)
    upsampler_model.eval()
    upsampler_diffusion = diffusion_from_config(DIFFUSION_CONFIGS['upsample'])

    print('downloading base checkpoint...')
    base_model.load_state_dict(load_checkpoint(base_name, device))

    print('downloading upsampler checkpoint...')
    upsampler_model.load_state_dict(load_checkpoint('upsample', device))

    pc_img_folder = f'images/pc_dalle/{base_name}/img'
    pc_npy_folder = f'images/pc_dalle/{base_name}/npy'
    os.makedirs(pc_img_folder, exist_ok=True)
    os.makedirs(pc_npy_folder, exist_ok=True)

    sampler = PointCloudSampler(
        device=device,
        models=[base_model, upsampler_model],
        diffusions=[base_diffusion, upsampler_diffusion],
        num_points=[1024, 4096 - 1024],
        aux_channels=['R', 'G', 'B'],
        guidance_scale=[3.0, 3.0],
    )
    
    for image_name in tqdm(os.listdir(dalle_folder)):
        img = Image.open(osp.join(dalle_folder, image_name))
        samples = None
        for x in tqdm(sampler.sample_batch_progressive(batch_size=1, model_kwargs=dict(images=[img]))):
            samples = x
            
        pc = sampler.output_to_point_clouds(samples)[0]
        fig = plot_point_cloud(pc, grid_size=3, dot_size=0.1, fixed_bounds=((-0.75, -0.75, -0.75), (0.75, 0.75, 0.75)))
        plt.savefig(osp.join(pc_img_folder, image_name))
        plt.close()
        
        pc_save = np.concatenate([pc.coords, pc.channels['R'][:,None], pc.channels['G'][:,None], pc.channels['B'][:,None]], axis=1)
        np.save(osp.join(pc_npy_folder, f'{osp.splitext(image_name)[0]}'), pc_save)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--pointe_name', type=str, choices=['base300M', 'base1B'])
    parser.add_argument('--cleaned', action='store_true')
    args = parser.parse_args()    
    
    main(args)
    
"""   
python clean_images.py \
  && python dalle2pointcloud.py --device cuda:2 --pointe_name base1B \
  && python dalle2pointcloud.py --device cuda:2 --pointe_name base1B --cleaned \
  && python dalle2pointcloud.py --device cuda:2 --pointe_name base300M \
  && python dalle2pointcloud.py --device cuda:2 --pointe_name base300M --cleaned
""" 