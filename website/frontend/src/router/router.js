import Vue from 'vue'
import VueRouter from 'vue-router'
import Home from '../views/Home.vue'

Vue.use(VueRouter);

// Lazy loading
const About = () => import('../components/About.vue');
const Prediction = () => import('../components/Prediction.vue');
const AMPs = () => import('../components/AMPs.vue');

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
        title:'Macrel',
      },
      redirect: '/prediction',

      children:[
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
    base: '/software/macrel/',

});

export default router
