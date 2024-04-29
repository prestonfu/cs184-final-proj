import os
import os.path as osp

import torch
from tqdm import tqdm
from torch import autocast
from diffusers import StableDiffusionPipeline

sd_folder = 'images/stable_diffusion'
os.makedirs(sd_folder, exist_ok=True)

pipe = StableDiffusionPipeline.from_pretrained("stabilityai/stable-diffusion-2-1").to("cuda:7")
pipe.enable_attention_slicing()

with open('stable_diffusion_prompts.txt') as f:
    sd_prompts = [x.strip() for x in f.readlines()]

for prompt in tqdm(sd_prompts):
    fname = osp.join(sd_folder, f"{prompt.lower().replace(' ', '_')}.png")
    if osp.exists(fname): continue
    image = pipe(prompt).images[0]
    image.save(fname)