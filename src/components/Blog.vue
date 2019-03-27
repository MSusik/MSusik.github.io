
<template>
  <div>
    <article>
      <transition name="fade" mode="out-in">
      <router-view></router-view>
      </transition>
    </article>
  </div>
</template>

<script>
export default {
  name: 'Blog'
}
</script>

<style lang="scss" src="../assets/Blog.scss"/>
