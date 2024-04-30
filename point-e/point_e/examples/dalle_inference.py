import os
import os.path as osp
from dotenv import load_dotenv
from openai import OpenAI
from tqdm import tqdm
from urllib.request import urlretrieve
from urllib.error import HTTPError

load_dotenv()
client = OpenAI()

dalle_folder = 'images/dalle'
os.makedirs(dalle_folder, exist_ok=True)

with open('stable_diffusion_prompts.txt') as f:
    sd_prompts = [x.strip() for x in f.readlines()]

for prompt in tqdm(sd_prompts):
# for _ in range(1):
    fname = osp.join(dalle_folder, f"{prompt.lower().replace(' ', '_')}.png")
    # fname = osp.join(dalle_folder, 'blah.png')
    if osp.exists(fname): continue
    response = client.images.generate(
        model="dall-e-3",
        prompt=f"With a white background: {prompt}",
        size="1024x1024",
        quality="standard",
        n=1,
    )
    url = response.data[0].url
    try:
        urlretrieve(url, fname)
    except HTTPError:
        print(f'Could not retrieve {prompt} from {url}')