// 配置路由相关信息
import Vue from 'vue'
import VueRouter from 'vue-router'
import Home from '../views/Home.vue'

// 1.通过Vue.use(插件)，安装插件
Vue.use(VueRouter);

// 懒加载
const About = () => import('../components/About.vue');
const Prediction = () => import('../components/Prediction.vue');
const AMPs = () => import('../components/AMPs.vue');

// 配置路径和组件之间的关系
const routes = [

    {
      path: '/',
      name: 'index',
      redirect: '/home',
    },

    {
      path: '/home',
      name: 'home',
      component: Home,
      meta:{
        title:'FACS',
      },
      redirect: '/prediction',

      children:[
        {
          path:'/prediction',
          name:'prediction',
          component:Prediction,
          meta:{
            title:'prediction',
          },
        },

        {
          path:'/prediction/amps',
          name:'amps',
          component:AMPs,
          meta:{
            title:'results',
          },
        },

        {
          path:'/about',
          name:'about',
          component:About,
          meta:{
            title:'about',
          },
        },

      ],

    },

];

// 2.创建路由对象
const router = new VueRouter({
  routes,
  mode: 'history'
});

// 3.将router对象传入到Vue实例中
export default router
