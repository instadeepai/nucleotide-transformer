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

"""Package logger.

Library rules: attach only a ``NullHandler`` and let the application decide
where records go. ``cli.py`` calls ``logging.basicConfig``, which is what makes
a whole-genome run report progress instead of looking hung for 45 minutes; a
caller using ``run()`` programmatically configures logging as they please.

Attaching a real handler here instead would bind ``sys.stderr`` at import time,
which duplicates every line once the CLI installs a root handler and writes to
a closed stream under pytest's capture.
"""

import logging

_ROOT = "nucleotide_transformer_v3"


def get_logger(name: str = _ROOT) -> logging.Logger:
    """Return the package logger."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.addHandler(logging.NullHandler())
    return logger
