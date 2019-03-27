<template>
<div>
<h2>About Upmath</h2>
<p>It works in browsers, except equations rendered <a href="//tex.s2cms.com/">on the server</a>. The editor stores your text in the browser to prevent the loss of your work in case of software or hardware failures.</p>
<p>I have designed and developed this lightweight editor and the service for converting LaTeX equations into svg-pictures to make publishing math texts on the web easy. I consider client-side rendering, the rival technique implemented in <a href="https://www.mathjax.org/">MathJax</a>, to be too limited and resource-consuming, especially on mobile devices.</p>
<p>The source code is <a href="https://github.com/parpalak/upmath.me">published on Github</a> under MIT license.</p>
<hr>
<p>Now you can erase this instruction and start writing your own scientific post. If you want to see the instruction again, open the editor in a private tab, in a different browser or download and clear your post and refresh the page.</p>
<p>Have a nice day :)</p>
</div>
</template>

<script>
export default {
  name: 'About',
  data () {
    return {
      name: this.$route.params.id
    }
  }
}
</script>
