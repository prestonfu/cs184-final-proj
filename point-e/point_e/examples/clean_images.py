import os
import os.path as osp
import subprocess
from PIL import Image
from tqdm import tqdm

def clean_images(dirs, postprocessed_folder='images/postprocessed'):
    assert isinstance(dirs, list)
    all_images = []
    for d in dirs:
        all_images.extend([(d, x) for x in os.listdir(d)])
    os.makedirs(postprocessed_folder, exist_ok=True)
    
    for folder, image_name in tqdm(all_images):
        tmp_path = osp.join(postprocessed_folder, 'tmp_' + image_name)
        final_path = osp.join(postprocessed_folder, image_name)
        if osp.exists(final_path): continue
        
        # Recrop
        image = Image.open(osp.join(folder, image_name))
        width, height = image.size
        square_size = int(1.2 * max(width, height))
        new_image = Image.new("RGB", (square_size, square_size), color="white")

        x = (square_size - width) // 2
        y = (square_size - height) // 2
        new_image.paste(image, (x, y))

        # Make background transparent
        resized_image = new_image.resize((256, 256), Image.LANCZOS)
        resized_image.save(tmp_path)
        result = subprocess.run(
            ['backgroundremover', '-i', f'{tmp_path}', '-o', f'{final_path}'],
            capture_output=True,
            text=True
        )
        if result.stderr:
            print(result.stderr)
        
        # Add back white background
        image = Image.open(final_path)
        background = Image.new("RGB", image.size, (255, 255, 255))
        background.paste(image, mask=image.split()[3])  # 3 is the alpha channel
        background.save(final_path)
        os.remove(tmp_path)

if __name__ == '__main__':
    clean_images(['images/dalle'], 'images/postprocessed_dalle')