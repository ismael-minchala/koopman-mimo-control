/* =====================================================
   KOOPMAN MIMO RESEARCH · UC-DEET
   main.js — canvas, scroll reveal, counters, nav
   ===================================================== */

// ── NAVBAR ──────────────────────────────────────────────
const navbar    = document.getElementById('navbar');
const navToggle = document.getElementById('navToggle');
const navLinks  = document.getElementById('navLinks');

window.addEventListener('scroll', () => {
  navbar.classList.toggle('scrolled', window.scrollY > 60);
  updateActiveLink();
}, { passive: true });

navToggle.addEventListener('click', () => {
  navLinks.classList.toggle('open');
  navToggle.classList.toggle('open');
});
navLinks.querySelectorAll('a').forEach(a => {
  a.addEventListener('click', () => {
    navLinks.classList.remove('open');
    navToggle.classList.remove('open');
  });
});

function updateActiveLink() {
  const sections  = document.querySelectorAll('section[id]');
  const anchors   = navLinks.querySelectorAll('a');
  let   current   = '';
  sections.forEach(s => {
    if (window.scrollY >= s.offsetTop - 130) current = s.id;
  });
  anchors.forEach(a => {
    const active = a.getAttribute('href') === `#${current}`;
    a.style.color = active ? 'var(--cyan)' : '';
  });
}

// ── HERO CANVAS — Koopman particle field ─────────────────
(function initCanvas() {
  const canvas = document.getElementById('heroCanvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  const COLORS = ['rgba(6,182,212,', 'rgba(129,140,248,', 'rgba(52,211,153,'];
  let W, H, particles = [];

  function resize() {
    W = canvas.width  = window.innerWidth;
    H = canvas.height = window.innerHeight;
  }

  class Particle {
    constructor() { this.init(); }
    init() {
      this.x  = Math.random() * W;
      this.y  = Math.random() * H;
      this.r  = Math.random() * 1.8 + .6;
      this.vx = (Math.random() - .5) * .6;
      this.vy = (Math.random() - .5) * .6;
      this.c  = COLORS[Math.floor(Math.random() * COLORS.length)];
      this.a  = Math.random() * .5 + .15;
    }
    update() {
      this.x += this.vx;
      this.y += this.vy;
      if (this.x < 0) this.x = W;
      if (this.x > W) this.x = 0;
      if (this.y < 0) this.y = H;
      if (this.y > H) this.y = 0;
    }
    draw() {
      ctx.beginPath();
      ctx.arc(this.x, this.y, this.r, 0, Math.PI * 2);
      ctx.fillStyle = this.c + this.a + ')';
      ctx.fill();
    }
  }

  function buildParticles() {
    const count = Math.min(Math.floor(W * H / 14000), 90);
    particles = Array.from({ length: count }, () => new Particle());
  }

  function drawGrid() {
    ctx.strokeStyle = 'rgba(6,182,212,0.035)';
    ctx.lineWidth   = 1;
    const g = 70;
    for (let x = 0; x < W; x += g) {
      ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke();
    }
    for (let y = 0; y < H; y += g) {
      ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke();
    }
  }

  function drawLinks() {
    const MAX = 160;
    for (let i = 0; i < particles.length; i++) {
      for (let j = i + 1; j < particles.length; j++) {
        const dx   = particles[i].x - particles[j].x;
        const dy   = particles[i].y - particles[j].y;
        const dist = Math.hypot(dx, dy);
        if (dist < MAX) {
          const alpha = (1 - dist / MAX) * .28;
          ctx.beginPath();
          ctx.moveTo(particles[i].x, particles[i].y);
          ctx.lineTo(particles[j].x, particles[j].y);
          ctx.strokeStyle = `rgba(6,182,212,${alpha})`;
          ctx.lineWidth   = .5;
          ctx.stroke();
        }
      }
    }
  }

  // Animated eigenfunction wave overlay
  let t = 0;
  function drawWaves() {
    ctx.save();
    ctx.globalAlpha = .04;
    for (let k = 0; k < 3; k++) {
      ctx.beginPath();
      const amp  = 40 + k * 15;
      const freq = .008 - k * .002;
      const off  = (t * (.3 + k * .1)) % (Math.PI * 2);
      ctx.moveTo(0, H / 2);
      for (let x = 0; x <= W; x += 4) {
        const y = H / 2 + Math.sin(x * freq + off) * amp * Math.cos(x * freq * .5 + off * .7);
        ctx.lineTo(x, y);
      }
      ctx.strokeStyle = k === 0 ? '#06b6d4' : k === 1 ? '#818cf8' : '#34d399';
      ctx.lineWidth   = 1.5;
      ctx.stroke();
    }
    ctx.restore();
  }

  function loop() {
    ctx.clearRect(0, 0, W, H);
    t += .018;
    drawGrid();
    drawWaves();
    drawLinks();
    particles.forEach(p => { p.update(); p.draw(); });
    requestAnimationFrame(loop);
  }

  resize();
  buildParticles();
  loop();

  window.addEventListener('resize', () => { resize(); buildParticles(); }, { passive: true });
})();

// ── COUNTER ANIMATION ────────────────────────────────────
function animateCounter(el) {
  if (el.dataset.counted) return;
  el.dataset.counted = '1';
  const target = parseInt(el.dataset.target, 10);
  const suffix = el.dataset.suffix || '';
  const dur    = 1600;
  const start  = performance.now();
  function tick(now) {
    const p = Math.min((now - start) / dur, 1);
    const e = 1 - Math.pow(1 - p, 3);           // ease-out cubic
    el.textContent = Math.round(e * target) + suffix;
    if (p < 1) requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}

// Trigger hero counters after short delay (always visible)
setTimeout(() => {
  document.querySelectorAll('.hstat-num[data-target]').forEach(animateCounter);
}, 600);

// ── SCROLL REVEAL + COUNTER TRIGGER ─────────────────────
const io = new IntersectionObserver((entries) => {
  entries.forEach(entry => {
    if (!entry.isIntersecting) return;
    entry.target.classList.add('visible');
    // counters inside this element
    entry.target.querySelectorAll('[data-target]').forEach(animateCounter);
    // if the element itself is a counter
    if (entry.target.hasAttribute('data-target')) animateCounter(entry.target);
    io.unobserve(entry.target);
  });
}, { threshold: 0.12, rootMargin: '0px 0px -40px 0px' });

// Observe all reveal elements
document.querySelectorAll('.reveal').forEach(el => io.observe(el));

// Observe result counters directly
document.querySelectorAll('.rs-num[data-target]').forEach(el => io.observe(el));

// ── STAGGERED CHILDREN ───────────────────────────────────
const staggerParents = [
  '.highlights', '.repo-grid', '.method-grid',
  '.pub-list', '.team-grid', '.rq-grid', '.obj-list'
];
const staggerIo = new IntersectionObserver((entries) => {
  entries.forEach(entry => {
    if (!entry.isIntersecting) return;
    const children = entry.target.querySelectorAll('.reveal:not(.visible)');
    children.forEach((child, i) => {
      setTimeout(() => {
        child.classList.add('visible');
        child.querySelectorAll('[data-target]').forEach(animateCounter);
      }, i * 90);
    });
    staggerIo.unobserve(entry.target);
  });
}, { threshold: 0.08 });

staggerParents.forEach(sel => {
  document.querySelectorAll(sel).forEach(el => staggerIo.observe(el));
});

// ── CLOSE MOBILE NAV ON OUTSIDE CLICK ───────────────────
document.addEventListener('click', e => {
  if (navLinks.classList.contains('open') &&
      !navLinks.contains(e.target) &&
      !navToggle.contains(e.target)) {
    navLinks.classList.remove('open');
    navToggle.classList.remove('open');
  }
});
