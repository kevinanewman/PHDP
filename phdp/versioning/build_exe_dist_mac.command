#! /bin/zsh

cd ..

# build executable

pyinstaller phdp.py \
    --name PHDP-0.1.2-mac-arm64.command \
    --paths .:common \
    --add-data report_templates:report_templates \
    --noconfirm \
    --onefile

# cleanup

mv *.spec versioning
rm -R __pycache__
rm -R build

cd versioning
