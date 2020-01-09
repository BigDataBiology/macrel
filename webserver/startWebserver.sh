#!/bin/bash

#================================================================
#	author:HiramHe
#	description:start nginx.run jar file.
#================================================================


echo ''
echo '************************************************************************************************************'
echo "	Before starting webserver,please make sure that you hava installed nginx and java environment."
echo "	This script needs sudo privilege,or it won't work well."
echo '************************************************************************************************************'
echo ''


workdir=$(cd $(dirname $0); pwd)
nginx_path="$HOME/software/nginx"
pipelineHome=$(cd ..; pwd)
logFile='website.log'
#HOST_IP=$(ifconfig -a|grep inet|grep -v 127.0.0.1|grep -v inet6|awk '{print $2}'|tr -d "addr:")
HOST_IP=$(curl -s icanhazip.com)


if [ -z $HOST_IP ]
then
	HOST_IP="your internet ip"
fi


echo '========================start nginx============================'
nginxPID=$( ps -ef | grep nginx | grep -v grep | grep 'master process' | awk '{print $2}' )
if [ ! -z $nginxPID ]
then
	echo 'nginx is already running.'
	echo 'no action.'
else
	if [ -x "$nginx_path/nginx/sbin/nginx" ]
	then
		echo 'nginx installed.'
		echo 'starting...'
	
		sudo $nginx_path/nginx/sbin/nginx
	
		if [ $? -eq 0 ]
		then
			echo 'start successfully.'
			echo 'you can visit the webisite in browser with public network ip of your server.'
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

apiServicePID=$( ps -ef | grep java | grep -v grep | grep 'FACS' | awk '{print $2}' )
if [ ! -z $apiServicePID ]
then
	echo 'webserver is already running.'
	echo 'no action.'
else
	jarName=$(ls | grep ".jar$")
	if [ -z $jarName ]
	then
		echo 'jar file not found.process stop.'
		echo 'exit!'
		exit 1
	else
		echo 'jar file found.'
		echo 'start websiteEnd service by daemon...'
		
		nohup java -jar $jarName --pipeline.home=$pipelineHome >$logFile &
		
		if [ $? -eq 0 ]
		then
			echo 'start successfully.'
		else
			echo 'failed to start,exit.'
			exit 1
		fi
	fi
fi

echo '==============================starting websiteEnd finish==================='
echo ''

