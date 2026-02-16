// Project filtering functionality
document.addEventListener('DOMContentLoaded', function() {
  const filterButtons = document.querySelectorAll('.filter-btn');
  const projectCards = document.querySelectorAll('.project-card[data-domain]');

  if (filterButtons.length === 0 || projectCards.length === 0) return;

  filterButtons.forEach(function(btn) {
    btn.addEventListener('click', function() {
      // Update active button
      filterButtons.forEach(function(b) { b.classList.remove('active'); });
      btn.classList.add('active');

      const filter = btn.getAttribute('data-filter');

      projectCards.forEach(function(card) {
        if (filter === 'all') {
          card.style.display = '';
        } else {
          var domains = card.getAttribute('data-domain') || '';
          if (domains.toLowerCase().indexOf(filter.toLowerCase()) !== -1) {
            card.style.display = '';
          } else {
            card.style.display = 'none';
          }
        }
      });
    });
  });
});
