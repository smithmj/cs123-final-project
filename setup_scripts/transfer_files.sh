#!/bin/bash

cd 
#scp -i seq_align.pem cs123-final-project/seq_alignment/*.py ec2-user@$1

scp -i seq_align.pem cs123-final-project/setup.sh ec2-user@$1

#echo "does this work"$1 "and this"
