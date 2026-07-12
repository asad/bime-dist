/*
 * BIME — BioInception Molecular Editor
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK
 * Algorithm Copyright (c) 2026 Syed Asad Rahman
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
// v2.4.5: apply the colour-blind-safe (CVD) palette preference as early as
// nav.js executes — BEFORE DOMContentLoaded and before the workbench renders
// any coloured pathway/trace content — so a returning CVD user never sees a
// flash of the default palette. Mirrors the dark-mode `data-theme` mechanism
// (localStorage key + an attribute on <html>); the two are orthogonal and
// compose. localStorage may throw in private-browsing / sandboxed iframes, so
// a read failure or a foreign value must not block page init.
(function applyStoredPalette() {
    var stored = null;
    try { stored = localStorage.getItem('bime-palette'); } catch (e) { stored = null; }
    if (stored === 'cvd') {
        document.documentElement.setAttribute('data-palette', 'cvd');
    }
})();
document.addEventListener('DOMContentLoaded', function() {
    var btn = document.querySelector('.nav-hamburger');
    var links = document.querySelector('.nav-links');

    // Hamburger menu toggle. Mirror visual state to aria-expanded so
    // screen-reader users hear the same open/closed state. WCAG 4.1.2.
    if (btn && links) {
        var setExpanded = function(open) {
            btn.setAttribute('aria-expanded', open ? 'true' : 'false');
        };
        btn.addEventListener('click', function(e) {
            e.stopPropagation();
            var willOpen = !links.classList.contains('open');
            links.classList.toggle('open');
            setExpanded(willOpen);
        });
        // Close menu when clicking a link
        links.querySelectorAll('a').forEach(function(a) {
            a.addEventListener('click', function() {
                links.classList.remove('open');
                setExpanded(false);
            });
        });
        // Close menu when clicking outside
        document.addEventListener('click', function(e) {
            if (!btn.contains(e.target) && !links.contains(e.target)) {
                if (links.classList.contains('open')) {
                    links.classList.remove('open');
                    setExpanded(false);
                }
            }
        });
        // Close menu on Escape key (keyboard users escape modals/menus)
        document.addEventListener('keydown', function(e) {
            if (e.key === 'Escape' && links.classList.contains('open')) {
                links.classList.remove('open');
                setExpanded(false);
                btn.focus();
            }
        });
    }

    // Dark mode toggle
    var themeBtn = document.querySelector('.theme-toggle');
    if (themeBtn) {
        // localStorage may throw in private-browsing or sandboxed iframes,
        // and a foreign value should not blow up theme init.
        var stored = null;
        try { stored = localStorage.getItem('bime-theme'); } catch (e) { stored = null; }
        if (stored === 'dark' || stored === 'light') {
            document.documentElement.setAttribute('data-theme', stored);
        }
        updateThemeIcon(themeBtn);

        themeBtn.addEventListener('click', function() {
            var current = document.documentElement.getAttribute('data-theme');
            var isDark = current === 'dark' ||
                (!current && window.matchMedia('(prefers-color-scheme: dark)').matches);
            var next = isDark ? 'light' : 'dark';
            document.documentElement.setAttribute('data-theme', next);
            try { localStorage.setItem('bime-theme', next); } catch (e) { /* quota / disabled */ }
            updateThemeIcon(themeBtn);
        });
    }

    // v2.4.5: colour-blind-safe (CVD) palette toggle. Mirrors the dark-mode
    // toggle above — same localStorage mechanism (key 'bime-palette'), an
    // attribute on <html> (data-palette="cvd"), aria-pressed + a dynamic
    // aria-label for WCAG 4.1.2, and keyboard accessibility (it is a native
    // <button>, so Enter/Space and the global focus ring come for free). The
    // CVD palette is a binary on/off and composes with dark mode. The stored
    // preference was already applied before DOMContentLoaded by
    // applyStoredPalette() above; here we just sync the button + wire clicks.
    var paletteBtn = document.querySelector('.palette-toggle');
    if (paletteBtn) {
        updatePaletteIcon(paletteBtn);
        paletteBtn.addEventListener('click', function() {
            var isCvd = document.documentElement.getAttribute('data-palette') === 'cvd';
            if (isCvd) {
                document.documentElement.removeAttribute('data-palette');
                try { localStorage.setItem('bime-palette', 'default'); } catch (e) { /* quota / disabled */ }
            } else {
                document.documentElement.setAttribute('data-palette', 'cvd');
                try { localStorage.setItem('bime-palette', 'cvd'); } catch (e) { /* quota / disabled */ }
            }
            updatePaletteIcon(paletteBtn);
        });
    }

    function updateThemeIcon(el) {
        var current = document.documentElement.getAttribute('data-theme');
        var isDark = current === 'dark' ||
            (!current && window.matchMedia('(prefers-color-scheme: dark)').matches);
        // Reflect the toggle state and label dynamically. WCAG 4.1.2.
        el.setAttribute('aria-pressed', isDark ? 'true' : 'false');
        el.setAttribute('aria-label', isDark ? 'Switch to light mode' : 'Switch to dark mode');
        el.innerHTML = isDark
            ? '<svg aria-hidden="true" focusable="false" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/></svg>'
            : '<svg aria-hidden="true" focusable="false" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';
    }

    // v2.4.5: reflect the CVD-palette toggle state + label. WCAG 4.1.2 — the
    // pressed state and an action-describing label both update so screen-reader
    // and keyboard users hear/see the same on/off the .palette-toggle styling
    // shows visually. The icon is the standard contrast glyph (a circle with a
    // half-filled side) — a visual hint, marked aria-hidden so it adds no noise.
    function updatePaletteIcon(el) {
        var isCvd = document.documentElement.getAttribute('data-palette') === 'cvd';
        el.setAttribute('aria-pressed', isCvd ? 'true' : 'false');
        el.setAttribute('aria-label', isCvd
            ? 'Disable colour-blind-safe palette'
            : 'Enable colour-blind-safe palette');
        el.innerHTML = '<svg aria-hidden="true" focusable="false" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="9"/><path d="M12 3a9 9 0 0 0 0 18z" fill="currentColor" stroke="none"/></svg>';
    }
});
