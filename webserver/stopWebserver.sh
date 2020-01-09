#!/bin/bash

#==================================================================
#	author:HiramHe
#	description:stop nginx.stop webserver
#==================================================================


echo ''
echo '************************************************************************************************************'
echo "	This script needs sudo privilege."
echo '************************************************************************************************************'
echo ''


workDir=$( cd $(dirname $0); pwd )
nginx_path="$HOME/software/nginx"
logFile='website.log'


echo '==============stop nginx===================='
nginxPID=$( ps -ef | grep nginx | grep -v grep | grep 'master process' | awk '{print $2}' )
if [ ! -z $nginxPID ]
then
	echo 'nginx is running.'
	echo 'stoping...'
	
	sudo $nginx_path/nginx/sbin/nginx -s stop
	
	echo 'done'
else
	echo "nginx isn't running."
fi
echo '==============finish==========================='
echo ''


echo '==================stop webservice===================='
apiServicePID=$( ps -ef | grep java | grep -v grep | grep 'FACS' | awk '{print $2}' )
if [ ! -z $apiServicePID ]
then
	echo 'webserver is running.'
	echo 'stoping...'
	
	kill -9 $apiServicePID
	rm -rf $workDir/$logFile
	
	echo 'done'
else
	echo "webserver isn't running."
fi
echo '=====================finish==============================='
echo ''

