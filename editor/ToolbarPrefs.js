/**
 * ToolbarPrefs.js — BIME v1.4.3 toolbar customization preferences.
 *
 * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 - see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Persists user-customized toolbar layout (which groups + buttons appear,
 * in which order) to localStorage. Read by MolEditor's `_buildToolbar` /
 * `_buildAtomBar` at editor instantiation, written by the Customize panel
 * when the user changes preferences.
 *
 * Storage shape (localStorage key: 'bime-toolbar-prefs-v1'):
 *
 *   {
 *     "version": 1,
 *     "groupOrder": ["draw", "bonds", "rings", "edit", "view", "reaction", "actions", "export"],
 *     "groups": {
 *       "draw":     { "hidden": false, "items": ["bond", "chain"] },
 *       "bonds":    { "hidden": false, "items": ["single", "double", "triple", "stereo"] },
 *       ...
 *     },
 *     "atomBar":   ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "H"]
 *   }
 *
 * Hidden groups have `hidden: true` and are dropped at build time. Hidden
 * items are simply omitted from the `items` array. Any group / item not
 * listed in prefs falls back to the canonical default order from
 * MolEditor's TOOLBAR_GROUPS / ATOM_BAR — so adding a new button in a
 * future BIME version automatically appears for users with stored prefs
 * (it lands at the end of its group; the user can move it via the
 * Customize panel afterward).
 *
 * API
 * ---
 *   ToolbarPrefs.load()                   -> prefs object | null
 *   ToolbarPrefs.save(prefs)              -> bool (true on success)
 *   ToolbarPrefs.reset()                  -> bool (clears localStorage entry)
 *   ToolbarPrefs.applyToGroups(groups, p) -> filtered+reordered groups copy
 *   ToolbarPrefs.applyToAtomBar(bar, p)   -> filtered+reordered atom-bar copy
 *   ToolbarPrefs.snapshot(groups, bar)    -> a fresh prefs object matching
 *                                            the canonical defaults (used by
 *                                            the Customize panel as the
 *                                            "current state" baseline before
 *                                            the user toggles anything).
 *   ToolbarPrefs.STORAGE_KEY              -> 'bime-toolbar-prefs-v1'
 *   ToolbarPrefs.version                  -> '1.6.0'
 */
(function(global) {
    'use strict';

    var STORAGE_KEY = 'bime-toolbar-prefs-v1';
    var SCHEMA_VERSION = 1;

    // -----------------------------------------------------------------------
    // load() / save() / reset()
    //
    // localStorage is unavailable in private-browsing modes and in Node tests
    // — every caller path tolerates a missing storage and falls back to the
    // canonical defaults.
    // -----------------------------------------------------------------------
    function getStorage() {
        try {
            return (global && global.localStorage) ? global.localStorage : null;
        } catch (e) { return null; }
    }

    function load() {
        var s = getStorage();
        if (!s) { return null; }
        var raw;
        try { raw = s.getItem(STORAGE_KEY); } catch (e) { return null; }
        if (!raw) { return null; }
        var parsed;
        try { parsed = JSON.parse(raw); } catch (e) { return null; }
        if (!validate(parsed)) { return null; }
        return parsed;
    }

    function save(prefs) {
        if (!validate(prefs)) { return false; }
        var s = getStorage();
        if (!s) { return false; }
        try { s.setItem(STORAGE_KEY, JSON.stringify(prefs)); return true; }
        catch (e) { return false; }
    }

    function reset() {
        var s = getStorage();
        if (!s) { return false; }
        try { s.removeItem(STORAGE_KEY); return true; }
        catch (e) { return false; }
    }

    // -----------------------------------------------------------------------
    // validate(prefs) — defensive parse of an arbitrary object so that a
    // corrupted localStorage entry never crashes the toolbar build. Returns
    // true if the object satisfies the v1 schema.
    // -----------------------------------------------------------------------
    function validate(p) {
        if (!p || typeof p !== 'object') { return false; }
        if (p.version !== SCHEMA_VERSION) { return false; }
        if (!Array.isArray(p.groupOrder)) { return false; }
        if (!p.groups || typeof p.groups !== 'object') { return false; }
        if (!Array.isArray(p.atomBar)) { return false; }
        // Group entries: each must have hidden:bool, items: string[], and
        // hiddenItems: string[] (explicit user-hidden buttons).
        for (var k in p.groups) {
            if (!p.groups.hasOwnProperty(k)) { continue; }
            var g = p.groups[k];
            if (!g || typeof g !== 'object') { return false; }
            if (typeof g.hidden !== 'boolean') { return false; }
            if (!Array.isArray(g.items)) { return false; }
            if (g.hiddenItems !== undefined && !Array.isArray(g.hiddenItems)) {
                return false;
            }
            for (var i = 0; i < g.items.length; i++) {
                if (typeof g.items[i] !== 'string') { return false; }
            }
            if (g.hiddenItems) {
                for (var ih = 0; ih < g.hiddenItems.length; ih++) {
                    if (typeof g.hiddenItems[ih] !== 'string') { return false; }
                }
            }
        }
        // groupOrder must be all strings.
        for (var j = 0; j < p.groupOrder.length; j++) {
            if (typeof p.groupOrder[j] !== 'string') { return false; }
        }
        // atomBar must be all strings.
        for (var m = 0; m < p.atomBar.length; m++) {
            if (typeof p.atomBar[m] !== 'string') { return false; }
        }
        return true;
    }

    // -----------------------------------------------------------------------
    // snapshot(canonicalGroups, canonicalAtomBar)
    //
    // Build a fresh prefs object matching the canonical defaults. The
    // Customize panel uses this as the starting point; once the user
    // toggles a checkbox, this snapshot is mutated and saved.
    // -----------------------------------------------------------------------
    function snapshot(canonicalGroups, canonicalAtomBar) {
        var prefs = {
            version: SCHEMA_VERSION,
            groupOrder: [],
            groups: {},
            atomBar: (canonicalAtomBar || []).slice()
        };
        for (var i = 0; i < (canonicalGroups || []).length; i++) {
            var g = canonicalGroups[i];
            prefs.groupOrder.push(g.id);
            prefs.groups[g.id] = {
                hidden: false,
                items: (g.items || []).map(function(it) { return it.id; }),
                hiddenItems: []
            };
        }
        return prefs;
    }

    // -----------------------------------------------------------------------
    // applyToGroups(canonicalGroups, prefs)
    //
    // Returns a copy of `canonicalGroups` filtered + reordered per `prefs`.
    // Any group not mentioned in prefs.groupOrder is appended in its
    // canonical order (so newly-added BIME versions don't lose buttons for
    // users with stored prefs). Hidden groups are dropped. Within a group,
    // items follow prefs.groups[id].items order; un-listed items append in
    // their canonical order. Items in prefs that no longer exist (e.g. a
    // button was removed in a newer BIME version) are silently dropped.
    // -----------------------------------------------------------------------
    function applyToGroups(canonicalGroups, prefs) {
        if (!Array.isArray(canonicalGroups)) { return []; }
        if (!prefs || !validate(prefs)) { return canonicalGroups.slice(); }

        // Index canonical groups by id.
        var byId = {};
        for (var i = 0; i < canonicalGroups.length; i++) {
            byId[canonicalGroups[i].id] = canonicalGroups[i];
        }
        var seen = {};
        var out = [];
        // First, the prefs-ordered groups.
        for (var k = 0; k < prefs.groupOrder.length; k++) {
            var gid = prefs.groupOrder[k];
            if (!byId[gid] || seen[gid]) { continue; }
            seen[gid] = true;
            var pg = prefs.groups[gid];
            if (pg && pg.hidden) { continue; }
            out.push(reorderItems(byId[gid], pg));
        }
        // Then, any canonical group not mentioned in prefs (added in newer
        // BIME version after the prefs were last saved).
        for (i = 0; i < canonicalGroups.length; i++) {
            var cg = canonicalGroups[i];
            if (seen[cg.id]) { continue; }
            // Group not in prefs: keep visible, items in canonical order.
            out.push(cg);
        }
        return out;
    }

    function reorderItems(canonicalGroup, prefsGroup) {
        if (!prefsGroup || !Array.isArray(prefsGroup.items)) { return canonicalGroup; }
        var byId = {};
        for (var i = 0; i < (canonicalGroup.items || []).length; i++) {
            byId[canonicalGroup.items[i].id] = canonicalGroup.items[i];
        }
        // Items the user has explicitly hidden via the Customize panel.
        var hiddenSet = {};
        if (Array.isArray(prefsGroup.hiddenItems)) {
            for (var hi = 0; hi < prefsGroup.hiddenItems.length; hi++) {
                hiddenSet[prefsGroup.hiddenItems[hi]] = true;
            }
        }
        var seen = {};
        var items = [];
        // 1. Items in user's preferred order (skipping any explicitly hidden).
        for (var k = 0; k < prefsGroup.items.length; k++) {
            var iid = prefsGroup.items[k];
            if (!byId[iid] || seen[iid] || hiddenSet[iid]) { continue; }
            seen[iid] = true;
            items.push(byId[iid]);
        }
        // 2. Forward-compat: any canonical item NOT in prefs.items AND NOT
        // in prefs.hiddenItems is treated as a newer-BIME-version addition
        // and appended visible. The user can hide it via the Customize panel
        // afterwards (which moves it to hiddenItems).
        for (i = 0; i < (canonicalGroup.items || []).length; i++) {
            var ci = canonicalGroup.items[i];
            if (seen[ci.id] || hiddenSet[ci.id]) { continue; }
            items.push(ci);
        }
        // Shallow-copy the group, with the new items list.
        var copy = {};
        for (var key in canonicalGroup) {
            if (canonicalGroup.hasOwnProperty(key)) { copy[key] = canonicalGroup[key]; }
        }
        copy.items = items;
        return copy;
    }

    // -----------------------------------------------------------------------
    // applyToAtomBar(canonicalAtomBar, prefs)
    //
    // Returns a filtered+reordered atom-bar array. Symbols not in prefs are
    // dropped (the user explicitly didn't pick them). If prefs is missing
    // or invalid, returns the canonical bar unchanged.
    // -----------------------------------------------------------------------
    function applyToAtomBar(canonicalAtomBar, prefs) {
        if (!Array.isArray(canonicalAtomBar)) { return []; }
        if (!prefs || !validate(prefs)) { return canonicalAtomBar.slice(); }
        if (!Array.isArray(prefs.atomBar)) { return canonicalAtomBar.slice(); }
        // Restrict to symbols that exist in the canonical bar (prevents
        // injecting arbitrary strings via tampered localStorage).
        var canonSet = {};
        for (var i = 0; i < canonicalAtomBar.length; i++) { canonSet[canonicalAtomBar[i]] = true; }
        var out = [];
        var seen = {};
        for (var j = 0; j < prefs.atomBar.length; j++) {
            var sym = prefs.atomBar[j];
            if (canonSet[sym] && !seen[sym]) {
                out.push(sym);
                seen[sym] = true;
            }
        }
        return out;
    }

    var ToolbarPrefs = {
        load: load,
        save: save,
        reset: reset,
        validate: validate,
        snapshot: snapshot,
        applyToGroups: applyToGroups,
        applyToAtomBar: applyToAtomBar,
        STORAGE_KEY: STORAGE_KEY,
        SCHEMA_VERSION: SCHEMA_VERSION,
        version: '1.8.1'
    };

    if (typeof module !== 'undefined' && module.exports) {
        module.exports = ToolbarPrefs;
    }
    if (global) { global.ToolbarPrefs = ToolbarPrefs; }
})(typeof window !== 'undefined' ? window : (typeof global !== 'undefined' ? global : this));
