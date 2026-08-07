pip install --upgrade ..
python generate_slim.py > test_script.slim
sed -i '' '25s|.*|    defineConstant("trees_file", "/var/folders/lr/bq1yt55n7z97zchjxn12pyr40000gp/T/stdpopsim_o2acl0mp/c30ec0.trees")\;|' test_script.slim
