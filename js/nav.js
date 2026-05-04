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
});
