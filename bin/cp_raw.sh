cd /blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/raw

cp -r /orange/fengyue/CUT_RUN_JR/NS-4295-FYue .
cp -r /orange/fengyue/CUT_RUN_JR/NS-4363-FYue-10B-235T75LT3-Lane8 .
find ./NS-4192-FYue-10B-232CVGLT3-Lane8/  -type f -exec md5sum {} + > my_checksums.md5
find ./NS-4205-FYue-10B-232CVGLT3-Lane8/  -type f -exec md5sum {} + >> my_checksums.md5
find /orange/fengyue/NS-4192-FYue-10B-232CVGLT3-Lane8/ -type f -exec md5sum {} + > feng_checksums.md5
find /orange/fengyue/CUT_RUN_YZ_YW/NS-4205-FYue-10B-232CVGLT3-Lane8/ -type f -exec md5sum {} + >> feng_checksums.md5

