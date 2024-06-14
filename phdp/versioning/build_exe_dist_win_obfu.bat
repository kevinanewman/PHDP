REM usage build_dist

cd ..\..

rmdir /S /Q obfu
pyarmor gen -O obfu -r phdp
move /Y obfu\pyarmor_runtime_000000 obfu\phdp\pyarmor_runtime_000000
copy /Y phdp\report_templates obfu\phdp\report_templates

cd phdp\versioning
