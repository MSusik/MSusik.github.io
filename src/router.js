import Vue from 'vue'
import Router from 'vue-router'
import Upmath from '@/components/articles/Upmath'
import About from '@/components/articles/About'

Vue.use(Router)

export default new Router({
  routes: [
    {
      path: '/',
      redirect: '/upmath'
    },
    {
      path: '/upmath',
      name: 'Upmath',
      component: Upmath
    },
    {
      path: '/about',
      name: 'About',
      component: About
    }
  ]
})
