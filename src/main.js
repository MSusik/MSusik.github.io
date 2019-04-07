// The Vue build version to load with the `import` command
// (runtime-only or standalone) has been set in webpack.base.conf with an alias.
import Vue from 'vue'
import Vuex from 'vuex'
import App from './App'
import './../node_modules/bulma/bulma.sass'
import { library } from '@fortawesome/fontawesome-svg-core'
import { faGithub } from '@fortawesome/free-brands-svg-icons'
import { FontAwesomeIcon } from '@fortawesome/vue-fontawesome'
import router from './router'
import hljs from 'highlight.js/lib/highlight'
import cpp from 'highlight.js/lib/languages/cpp'
import 'highlight.js/styles/pojoaque.css'

library.add(faGithub)
hljs.registerLanguage('cpp', cpp)
hljs.initHighlightingOnLoad()

Vue.component('font-awesome-icon', FontAwesomeIcon)
Vue.config.productionTip = false
Vue.use(Vuex)

Object.defineProperty(Vue.prototype, '$hljs', { value: hljs })

const store = new Vuex.Store({
  state: {
    base: 'https://raw.githubusercontent.com/MSusik/msusik.github.io/development/src/assets/backs/',
    image: '2.png',
    images: ['1.png', '2.png', '3.png', '4.png', '5.png', '6.png'],
    backgroundInTransition: false
  },
  getters: {
    base: state => state.base,
    images: state => state.images,
    image: state => state.image,
    backgroundInTransition: state => state.backgroundInTransition
  },
  mutations: {
    cacheImages: function (state) {
      state.images.forEach(function (img) {
        new Image().src = state.base + img
        // caches images, avoiding white flash between background replacements
      })
    },
    changeImage (state, { image }) {
      state.image = image
    },
    changeBackgroundInTransition (state, { value }) {
      state.backgroundInTransition = value
    }
  }
})

/* eslint-disable no-new */
new Vue({
  el: '#app',
  store,
  router,
  components: { App },
  template: '<App/>'
})
