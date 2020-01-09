#!/bin/bash

#================================================================
#	author:HiramHe
#	description:start nginx.run jar file.
#================================================================

echo ''


workdir=$(cd $(dirname $0); pwd)
nginx_path="$HOME/software/nginx"
pipelineHome=$(cd ..; pwd)
logFile='website.log'
#HOST_IP=$(ifconfig -a|grep inet|grep -v 127.0.0.1|grep -v inet6|awk '{print $2}'|tr -d "addr:")
HOST_IP=$(curl -s icanhazip.com)
enableCustomizedSetting="false"


usage ()
{
	echo '************************************************************************************************************'
	echo "	Before starting webserver,please make sure that you hava installed nginx and java environment."
	echo "	This script needs sudo privilege,or it won't work well."
	echo '************************************************************************************************************'
	echo ''
	
	echo "
	usage: bash startWebserver.sh [-c/--customized true/false]
	options:
	-c	If you want to start webserver with your customized setting in application.properties file, your option should be:true.
		Default:false
	"
}


while [[ $# -gt 0 ]]
do
	case $1 in
		-h|--help)
            usage
            exit 0
        ;;
		-c|--customized)
			enableCustomizedSetting=${2}
			if [ "$enableCustomizedSetting"x == "true"x ]
			then
				echo 'webserver will be started with customized setting.'
			elif [ "$enableCustomizedSetting"x == "false"x ]
			then
				echo 'webserver will be started with default setting.'
			else
				echo "the parameter '-c' didn't get correct value,please check it first."
				echo 'exit!'
				exit 1
			fi
		;;
	esac
	shift
done


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
		
		if [ $enableCustomizedSetting == "true" ]
		then
			echo 'starting with customized setting...'
			nohup java -jar $jarName >$logFile &
		else
			echo 'starting with default setting...'
			nohup java -jar $jarName --pipeline.home=$pipelineHome >$logFile &
		fi
		
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

