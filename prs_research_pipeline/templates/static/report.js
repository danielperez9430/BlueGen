
function toggleSection(id, header) {
    const body = document.getElementById(id);
    const arrow = document.getElementById(id + '_arrow');
    if (body.style.display === 'none') {
        body.style.display = 'block';
        arrow.textContent = '▼';
    } else {
        body.style.display = 'none';
        arrow.textContent = '▶';
    }
}
function expandAll() {
    document.querySelectorAll('.section-body').forEach(b => b.style.display = 'block');
    document.querySelectorAll('.section-arrow').forEach(a => a.textContent = '▼');
}
function collapseAll() {
    document.querySelectorAll('.section-body').forEach(b => b.style.display = 'none');
    document.querySelectorAll('.section-arrow').forEach(a => a.textContent = '▶');
}
