#!/bin/bash --login
set +euo pipefail
conda activate vsd-1.0
set -euo pipefail

echo "CMD $@"
exec $@
