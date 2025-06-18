#!/bin/bash
# convert_images_to_webp.sh

for img in *.jpg *.jpeg *.png; do
  [ -e "$img" ] || continue  # skip if no files match
  base="${img%.*}"
  cwebp -q 80 "$img" -o "${base}.webp"
done
