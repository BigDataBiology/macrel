#!/bin/bash

#================================================================
#	author:HiramHe
#	description:install node,nginx,java,maven.compile source code of vue,java.
#================================================================


workdir=$(cd $(dirname $0); pwd)
node_path="$HOME/software/node"
nginx_path="$HOME/software/nginx"
java_path="$HOME/software/java"
mvn_path="$HOME/software/maven"
echo ''


echo "===============================install npm======================================"
echo "check npm"
type npm >/dev/null 2>&1
if [ $? -eq 0 ]
then
  echo "npm installed"
else
  echo "npm not installed,we will install it"
  mkdir -p $node_path
  cd $node_path
  
  wget -O $node_path/node-v9.9.0-linux-x64.tar.xz https://nodejs.org/dist/v9.9.0/node-v9.9.0-linux-x64.tar.xz
  
  tar -xvf $node_path/node-v9.9.0-linux-x64.tar.xz -C $node_path >/dev/null 2>&1
  sudo ln -s $node_path/node-v9.9.0-linux-x64/bin/node /usr/bin/node
  sudo ln -s $node_path/node-v9.9.0-linux-x64/bin/npm /usr/bin/npm
  
  rm -rf $node_path/node-v9.9.0-linux-x64.tar.xz
  
fi
cd $workdir
echo "===========================finish to install npm====================================="
echo ''


echo '==================================install nginx======================================'
mkdir -p $nginx_path
cd $nginx_path

wget -O $nginx_path/openssl-1.0.2s.tar.gz https://www.openssl.org/source/openssl-1.0.2s.tar.gz
wget -O $nginx_path/pcre-8.43.tar.gz https://ftp.pcre.org/pub/pcre/pcre-8.43.tar.gz
wget -O $nginx_path/zlib-1.2.11.tar.gz https://zlib.net/zlib-1.2.11.tar.gz
wget -O $nginx_path/nginx-1.17.1.tar.gz http://nginx.org/download/nginx-1.17.1.tar.gz

echo "unzip nginx archive"
ls *.tar.gz | xargs -n1 tar xzvf >/dev/null 2>&1

cd $nginx_path/nginx-1.17.1
echo 'configure nginx'
./configure \
--with-openssl=../openssl-1.0.2s \
--with-pcre=../pcre-8.43 \
--with-zlib=../zlib-1.2.11 \
--prefix=$nginx_path/nginx \
--with-http_ssl_module \
--with-http_v2_module \
>/dev/null 2>&1

echo 'compile nginx'
make >/dev/null 2>&1
make install >/dev/null 2>&1
cd $nginx_path

echo 'remove package file'
ls *.tar.gz | xargs rm -rf
rm -rf nginx-1.17.1
rm -rf openssl-1.0.2s
rm -rf pcre-8.43
rm -rf zlib-1.2.11

echo "install nginx successfully"
cd $workdir
echo '===================================finish to install nginx=========================='
echo ''


echo '===================================install java==============================='
echo 'check java'
type java >/dev/null 2>&1
if [ $? -eq 0 ]
then
	echo "java was installed!"
else
	echo "java not installed!we are going to install it!"
	mkdir -p $java_path
	cd $java_path
	
	wget -O jdk-8u131-linux-x64.tar.gz \
	--no-check-certificate --no-cookies \
	--header "Cookie: oraclelicense=accept-securebackup-cookie" \
	http://download.oracle.com/otn-pub/java/jdk/8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz
	
	jdkName=$( ls | grep jdk-*.gz )
	if [ -f $jdkName ]
	then
		echo "upzip jdk-*-linux-x64.tar.gz"
		tar -zxvf $jdkName -C $java_path >/dev/null 2>&1
		
		echo "configure java environment."
		echo "#java jdk" >> ~/.bashrc
		echo "export JAVA_HOME=$java_path/jdk1.8.0_131" >> ~/.bashrc
		echo "export JRE_HOME=\$JAVA_HOME/jre" >> ~/.bashrc
		echo "export CLASSPATH=.:\$JAVA_HOME/lib:\$JRE_HOME/lib" >> ~/.bashrc
		echo "export PATH=\$JAVA_HOME/bin:\$PATH" >> ~/.bashrc
		
		echo "installing java finish."
		echo "you can use this command to activate it,or you can not use java."
		echo "source ~/.bashrc"
		
		rm -rf $jdkName
		
	else
		echo 'jdk file not exists,exit.'
		exit 1
	fi
fi
cd $workdir
echo '==========================finish to install java=============================='
echo ''


echo '============================install maven================================'
echo "check maven"
type mvn >/dev/null 2>&1
if [ $? -eq 0 ]
then
	echo "maven installed."
else
	echo "maven not installed.we are going to install it."
	mkdir -p $mvn_path
	cd $mvn_path
	
	echo "downloading maven package."
	wget -O apache-maven-3.3.9-bin.tar.gz  \
	--no-cookies --no-check-certificate \
	--header "Cookie: oraclelicense=accept-securebackup-cookie" \
	"http://mirrors.hust.edu.cn/apache/maven/maven-3/3.3.9/binaries/apache-maven-3.3.9-bin.tar.gz"
	
	mvnfile=$( ls | grep apache*maven-*.gz	)
	if [ -f $mvnfile ];
	then
		echo 'unzip maven archive.'
		tar zvxf $mvnfile -C $mvn_path --strip-components 1 >/dev/null 2>&1
		echo "unzip maven successfully"
		
		echo "configure environment variables"
		mv ~/.bashrc ~/.bashrc.backup.mvn
		cat ~/.bashrc.backup.mvn >> ~/.bashrc
		echo "# maven" >> ~/.bashrc
		echo "export MAVEN_HOME=$mvn_path" >> ~/.bashrc
		echo "export PATH=\$MAVEN_HOME/bin:\$PATH" >> ~/.bashrc
		echo "configure finish."
		
		echo "you can use this command to activate it,or you can not use mvn."
		echo "source ~/.bashrc"
		
		rm -rf $mvnfile
		
	else
		echo "no mvn archive.stop installation."
		exit 1
	fi	
fi
cd $workdir
echo '=========================finish to install maven===================================='
echo ''


echo 'install webserver successfully.'
echo "please use this command to activate environment:"
echo "source ~/.bashrc"
source ~/.bashrc
echo '===================================thanks for installing webserver==========================================='
