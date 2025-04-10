# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import os
from typing import Any, Dict, Tuple

import haiku as hk
import joblib
from huggingface_hub import hf_hub_download

from nucleotide_transformer.sCellTransformer.model import sCTConfig

ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"


def _get_dir(model_name: str) -> str:
    """
    Get directory to save files on user machine.
    """
    return os.path.expanduser(
        os.path.join(os.getenv(ENV_XDG_CACHE_HOME, DEFAULT_CACHE_DIR), model_name)
    )


def download_ckpt() -> Tuple[hk.Params, Any]:
    """
    Download checkpoint on kao datacenter.

    Args:
        model_name: Name of the model.

    Returns:
        Model parameters.
        Model state


    """

    save_dir = os.path.join(_get_dir("sCellTransformer"), "sCellTransformer")

    repo_id = f"InstaDeepAI/sCellTransformer"

    # Download parameters
    print("Downloading model's weights...")
    params = joblib.load(
        hf_hub_download(
            repo_id=repo_id,
            filename="jax_params/params.joblib",
            cache_dir=save_dir,
        )
    )

    config_path = hf_hub_download(
        repo_id=repo_id,
        filename="jax_params/config.json",
    )
    with open(config_path, "r") as f:
        config_dict = json.load(f)
    config = sCTConfig(**config_dict)

    return params, config
