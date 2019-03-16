
<template>
  <div>
    <div id="output"></div>

    <div v-for="n in 40" :key="n">
      {{ n }}
    </div>
  </div>
</template>

<script>
export default {
  name: 'Blog',
  methods: {
    displayArticle: function() {}
  },
  beforeMount () {
    this.displayArticle()
  },
}
</script>

<style lang="scss" src="../assets/Blog.scss"/>
