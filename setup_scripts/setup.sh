#!/bin/bash

sudo mkdir /mnt/align/py_code
sudo mount /dev/sdf /mnt/align

sudo mv *.py /mnt/align/pycode 

sudo yum install numpy 

sudo yum -y groupinstall "Development Tools"
sudo yum -y install fuse fuse-devel autoconf automake curl-devel libxml2-devel openssl-devel mailcap
wget https://github.com/s3fs-fuse/s3fs-fuse/archive/v1.78.tar.gz
tar xzf v1.78.tar.gz 
cd s3fs-fuse-1.78
./autogen.sh 
./configure
make
sudo make install
cd
echo AKIAJURJRV6YYMA7M7PA:U6rGhy2CanEhYR2oGCKlOfbt+sggifEOFg0+8u6d > .passwd-s3fs
chmod 600 .passwd-s3fs
mkdir new_data
s3fs cs123-smithmj-j3martinez-reads new_data

# to scp from vitual box to em2:

sudo chmod -R a+rwx /mnt/temp

# to install numpy on aws

sudo yum -y install gcc-c++ python27-devel atlas-sse3-devel lapack-devel
wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.2.tar.gz
tar xzf virtualenv-1.11.2.tar.gz 
python27 virtualenv-1.11.2/virtualenv.py sk-learn
. sk-learn/bin/activate
pip install numpy




