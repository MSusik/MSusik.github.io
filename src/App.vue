<template>
  <div class="wrap" id="app">
    <div class="left">
      <app-panel></app-panel>
    </div>
    <div class="right">
      <app-blog></app-blog>
    </div>
  </div>
</template>

<script>
import Panel from './components/Panel'
import Blog from './components/Blog'

export default {
  name: 'App',
  components: {
    'app-panel': Panel,
    'app-blog': Blog
  }
}
</script>

<style lang="scss" src="./assets/App.scss"/>
