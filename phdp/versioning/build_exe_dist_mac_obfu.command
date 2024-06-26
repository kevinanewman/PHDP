#! /bin/zsh

cd ..

# build executable and generate .spec file for pyarmor

pyinstaller phdp.py \
    --name PHDP-0.1.4-mac-arm64.command \
    --paths .:common \
    --add-data report_templates:report_templates \
    --noconfirm \
    --onefile

# generate obfuscated executable

pyarmor gen --pack PHDP-0.1.4-mac-arm64.command.spec phdp.py

# cleanup

mv *.spec versioning
rm -R .pyarmor
rm -R build

cd versioning
