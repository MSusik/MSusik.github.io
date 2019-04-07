<template>

  <main>
    <div class="panel1"/>
    <div class="panel2 content">

      <section>
      <div class="about">

        <div class="avatar"><img class="avatar" src="https://avatars1.githubusercontent.com/u/6649290?s=460&v=4" alt="avatar"></div>

        <h3>Mateusz Susik</h3>
        <h6>Machine learning developer</h6>

        <span class="icon iconlink">
          <a href="http://www.github.com/msusik">
            <i class="fab fa-github" href="http://www.github.com/msusik"></i>
          </a>
        </span>

        <span class="icon iconlink">
          <a href="http://www.linkedin.com/in/mateuszsusik">
            <i class="fab fa-linkedin"></i>
          </a>
        </span>

        <span class="icon iconlink">
          <a href="https://scholar.google.pl/citations?user=qAkjUM4AAAAJ&hl=en">
            <i class="fas fa-flask"></i>
          </a>
        </span>

        <p/>

        <span>Welcome to my blog on machine learning from various perspectives. I will cover here machine learning, cloud infrastructure and bioinformatics.</span>

      </div>
    </section>
  </div>
  <div class="panel3">
    <ul><li>
    <h6>Articles</h6>
    </li></ul>
    <ul>
      <li>
        <router-link :to="{ name: 'CudaPart1' }"><span v-on:click="changeRandomlyBackground">Implementing Smith-Waterman algorithm on CUDA. Part 1: C++ implementation</span></router-link>
      </li>
    </ul>
  </div>

  <div class="panel4">
     <small>Design and implementation - Mateusz Susik© 2019</small>
  </div>
  </main>

</template>

<script>

export default {
  name: 'Panel',
  methods: {
    changeRandomlyBackground: function () {
      if (!this.$store.getters.backgroundInTransition) {
        this.$store.commit('changeBackgroundInTransition', {value: true})
        var randImage = this.$store.getters.image
        while (randImage === this.$store.getters.image) {
          randImage = this.$store.getters.images[Math.floor(Math.random() * this.$store.getters.images.length)]
        }
        document.getElementsByClassName('bg-left')[0].style.background = 'url(\'' + this.$store.getters.base + randImage + '\')'
        this.$store.commit('changeImage', {image: randImage})
        setTimeout(function () {
          this.$store.commit('changeBackgroundInTransition', {value: false})
        }.bind(this), 1000)
      }
    }
  },
  created: function () {
    this.interval = setInterval(() => this.$store.commit('cacheImages'), 1000)
    this.$store.commit('cacheImages')
  }
}
</script>

<style lang="scss" src="../assets/Panel.scss"/>
