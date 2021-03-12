# UCAR|NCAR Icomoon icons
This repo is used as an npm package for Icomoon icons used in the custom Koru theme.

## Update Icomoon icons
To update the IcoMoon icon set (add, remove, etc):
 1. Go to the [IcoMoon online app](https://icomoon.io/app/#/select) (no need to login)
 2. Click "Import Icons" and choose "selection.json" from /source/fonts/icomoon/selection.json
 3. In the new set that is created, you can add, rearrange, deselect icons using the hamburger menu on the top right of the set
 4. For instance, to add a new icon, click the hamburger menu, then "Import to Set"
 5. Once the set is to your liking, click "Generate Font" in the lower right
 6. If the resulting page looks good, click "Download" on the same tab
 7. Unzip the download and overwrite the files in this repository

## Publish the icomoon package to Github Registry
Commit all new changes to Github. This will trigger a CircleCI build that will update the patch version, commit and push the changes to github (skipping the ci build in the process), packaging the code and publishing it to Github's registry.

To perform a minor update, run `npm version minor`.

To perform a major update, run `npm version major`.

After you've updated the minor or major version, commit and push the changes to Github.
