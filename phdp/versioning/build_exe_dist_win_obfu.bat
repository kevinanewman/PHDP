
cd ..

REM build executable and generate .spec file for pyarmor

pyinstaller phdp.py ^
    --name PHDP-0.1.5-win-obfu ^
    --paths .;common ^
    --add-data report_templates;report_templates ^
    --noconfirm ^
    --onefile

REM generate obfuscated executable

pyarmor gen --pack PHDP-0.1.5-win-obfu.spec phdp.py

REM cleanup

move /Y  *.spec versioning
rmdir /S /Q __pycache__
rmdir /S /Q build

cd versioning