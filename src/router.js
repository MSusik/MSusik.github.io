import Vue from 'vue'
import Router from 'vue-router'
import Upmath from '@/components/articles/Upmath'
import CudaPart1 from '@/components/articles/CudaPart1'

Vue.use(Router)

export default new Router({
  routes: [
    {
      path: '/',
      redirect: '/implementing_smith_waterman_on_CUDA_part_1'
    },
    {
      path: '/upmath',
      name: 'Upmath',
      component: Upmath
    },
    {
      path: '/implementing_smith_waterman_on_CUDA_part_1',
      name: 'CudaPart1',
      component: CudaPart1
    }
  ]
})
