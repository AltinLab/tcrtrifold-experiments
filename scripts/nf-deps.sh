#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <install|update>"
  exit 1
fi

action="$1"
if [[ "$action" != "install" && "$action" != "update" ]]; then
  echo "Error: Invalid action '$action'. Must be 'install' or 'update'."
  exit 1
fi

preview_flag=""
if [ "$action" = update ]; then
  preview_flag="--no-preview"
fi

nf-core subworkflows \
  --git-remote git@github.com:ljwoods2/af3-nf-tools.git \
  $action af3 --dir workflows $preview_flag

nf-core modules \
  --git-remote git@github.com:ljwoods2/af3-nf-tools.git \
  $action af3 --dir workflows $preview_flag