#! /bin/zsh

cd ../..

rm -R obfu
pyarmor gen -O obfu -r phdp
mv obfu/pyarmor_runtime_000000 obfu/phdp/pyarmor_runtime_000000
cp -R phdp/report_templates obfu/phdp/report_templates

cd obfu/phdp

# build executable

pyinstaller phdp.py \
    --name PHDP-0.1.4-mac-arm64.command \
    --paths .:common \
    --add-data report_templates:report_templates \
    --noconfirm \
    --onefile

# cleanup

mv *.spec versioning
rm -R __pycache__
rm -R build

cd versioning
