# 设置基础镜像
FROM nginx
# 定义作者
MAINTAINER HiramHe <hiram.he@qq.com>
# 将 dist 文件夹中的内容 复制到镜像的 /usr/share/nginx/html/ 自定义目录下，然后在nginx里面配好即可
COPY dist/ /usr/share/nginx/html/
COPY nginx.conf /etc/nginx/nginx.conf
RUN echo 'done!'
