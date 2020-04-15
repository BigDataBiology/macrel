import Vue from 'vue'
import VueRouter from 'vue-router'
import Home from '../views/Home.vue'
import Choice from "../components/Choice";

Vue.use(VueRouter);

// Lazy loading
const About = () => import('../components/About.vue');
const Prediction = () => import('../components/Prediction.vue');
const AMPs = () => import('../components/AMPs.vue');

const routes = [
    {
      //Equivalent to the root path at where app running.
      path: '/',
      name: 'index',
      redirect: '/home',
    },

    {
      path: '/home',
      name: 'home',
      component: Home,
      meta:{
        title:'Macrel',
      },
      redirect: '/choice',

      children:[
        {
          path:'/choice',
          name:'choice',
          component:Choice,
          meta:{
            title:'choose',
          },
        },

        {
          path:'/prediction',
          name:'prediction',
          component:Prediction,
          meta:{
            title:'AMP prediction',
          },
        },

        {
          path:'/prediction/amps',
          name:'amps',
          component:AMPs,
          meta:{
            title:'Results',
          },
        },

        {
          path:'/about',
          name:'about',
          component:About,
          meta:{
            title:'About',
          },
        },

      ],

    },

];

const router = new VueRouter({
    routes,
    mode: 'history',
    //equals to the value of publicPath in vue.config.js file.
    //So the router will jump correctly according to the root path at where app running.
    base: '/software/macrel/',

});

export default router
