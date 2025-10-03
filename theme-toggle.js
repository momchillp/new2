const toggle = document.getElementById('themeToggle');
const label = document.getElementById('theme-label');
const prefersDarkScheme = window.matchMedia("(prefers-color-scheme: dark)");

function applyTheme(theme) {
  if (theme === "dark") {
    document.documentElement.setAttribute("data-theme", "dark");
    label.textContent = "🌙 Dark";
    toggle.checked = true;
  } else {
    document.documentElement.setAttribute("data-theme", "light");
    label.textContent = "🌞 Light";
    toggle.checked = false;
  }
}

const storedTheme = localStorage.getItem("theme");
if (storedTheme) {
  applyTheme(storedTheme);
} else {
  applyTheme(prefersDarkScheme.matches ? "dark" : "light");
}

toggle.addEventListener("change", () => {
  const theme = toggle.checked ? "dark" : "light";
  applyTheme(theme);
  localStorage.setItem("theme", theme);
});
