import Vue from 'vue'
import App from './App.vue'
import router from './router/router'
import './plugins/element.js'
import axios from 'axios'
// 导入全局样式
import('./assets/css/global.css');

router.beforeEach( (to,from,next) =>{
  if (to.meta.title){
    /* 路由发生变化修改页面title */
    document.title=to.meta.title;
  }
  next();
} );

Vue.config.productionTip = false;

// 配置请求的根路径
axios.defaults.baseURL= 'http://39.106.68.204:8080/';
// axios.defaults.baseURL='http://localhost:8080';
// 定义成全局变量，所有vue组件都可以用
Vue.prototype.$http = axios;

new Vue({
  // el: '#app',
  router,
  axios,
  render: h => h(App)
}).$mount('#app');
