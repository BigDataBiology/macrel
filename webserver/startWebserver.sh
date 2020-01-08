#!/bin/bash

#================================================================
#	author:HiramHe
#	description:start nginx.run jar file.
#================================================================

echo '******************************************************************************************************'
echo "Before starting webserver,please make sure that you hava installed nginx and java environment."
echo "This script needs sudo privilege,or it won't work well."
echo '******************************************************************************************************'

workdir=$(cd $(dirname $0); pwd)
nginx_path="$HOME/software/nginx"
pipelineHome=$(cd ..; pwd)
HOST_IP=$(ifconfig -a|grep inet|grep -v 127.0.0.1|grep -v inet6|awk '{print $2}'|tr -d "addr:")


echo '========================start nginx============================'
if [ -d "$nginx_path/nginx" ]
then
	echo 'nginx installed.'
	echo 'start...'
	
	sudo $nginx_path/nginx/sbin/nginx
	
	if [ $? -eq 0 ]
	then
		echo 'start successfully.'
		echo 'you can visit the webisite in browser by public network ip of your server.'
		echo 'such as:'
		echo "http://$HOST_IP:80"
	else
		echo 'fail to start nginx.maybe you should reinstall it.exit!'
		exit 1
	fi
else
	echo 'nginx not installed.please install it first.exit!'
	exit 1
fi
echo '==============================starting nginx finish==================='
echo ''


echo '==============================start websiteEnd======================'
type java >/dev/null 2>&1
if [ $? -ne 0 ]
then
	echo 'java not installed.please install java first.'
	echo 'exit!'
	exit 1
fi

jarName=$(ls | grep ".jar$")
if [ -z $jarName ]
then
	echo 'jar file not found.process stop.'
	echo 'exit!'
	exit 1
else
	echo 'jar file found.'
	echo 'start websiteEnd service by daemon...'
	
	nohup java -jar $jarName --pipeline.home=$pipelineHome >website.log &
	
	if [ $? -eq 0 ]
	then
		echo 'start successfully.'
	else
		echo 'failed to start,exit.'
		exit 1
	fi
fi
echo '==============================starting websiteEnd finish==================='
echo ''

