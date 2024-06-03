REM usage build_dist

cd ..

REM build executable

pyinstaller phdp.py ^
    --name PHDP-0.1.1-win ^
    --paths .;common ^
    --add-data report_templates;report_templates ^
    --noconfirm ^
    --onefile

REM cleanup

move /Y  *.spec versioning
rmdir /S /Q __pycache__
rmdir /S /Q build

cd versioning
