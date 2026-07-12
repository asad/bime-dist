/**
 * History.js — Undo/redo stack for molecule editor
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 */
(function(global) {
    'use strict';

    var MAX_HISTORY = 100;

    // Ring-buffer-backed stack: avoids O(n) Array.shift() when the cap is hit
    // on every push.  Provides push / pop / length / clear semantics.
    function RingStack(cap) {
        this._buf = new Array(cap);
        this._cap = cap;
        this._head = 0;  // index of bottom-most element
        this._len = 0;   // number of elements currently stored
    }
    RingStack.prototype.push = function(item) {
        if (this._len < this._cap) {
            this._buf[(this._head + this._len) % this._cap] = item;
            this._len++;
        } else {
            // Overwrite the oldest entry; advance head
            this._buf[this._head] = item;
            this._head = (this._head + 1) % this._cap;
        }
    };
    RingStack.prototype.pop = function() {
        if (this._len === 0) return undefined;
        this._len--;
        var idx = (this._head + this._len) % this._cap;
        var item = this._buf[idx];
        this._buf[idx] = undefined; // release reference
        return item;
    };
    RingStack.prototype.clear = function() {
        // Drop references so the GC can reclaim large molecule snapshots
        for (var i = 0; i < this._cap; i++) this._buf[i] = undefined;
        this._head = 0;
        this._len = 0;
    };
    Object.defineProperty(RingStack.prototype, 'length', {
        get: function() { return this._len; }
    });

    function History() {
        this._undoStack = new RingStack(MAX_HISTORY);
        this._redoStack = new RingStack(MAX_HISTORY);
    }

    History.prototype.push = function(molecule) {
        this._undoStack.push(molecule.toJSON());
        // Pushing a new state invalidates the redo stack
        this._redoStack.clear();
    };

    History.prototype.undo = function(currentMolecule) {
        if (this._undoStack.length === 0) return null;
        this._redoStack.push(currentMolecule.toJSON());
        return this._undoStack.pop();
    };

    History.prototype.redo = function(currentMolecule) {
        if (this._redoStack.length === 0) return null;
        this._undoStack.push(currentMolecule.toJSON());
        return this._redoStack.pop();
    };

    History.prototype.canUndo = function() { return this._undoStack.length > 0; };
    History.prototype.canRedo = function() { return this._redoStack.length > 0; };

    History.prototype.clear = function() {
        this._undoStack.clear();
        this._redoStack.clear();
    };

    global.EditorHistory = History;

})(typeof window !== 'undefined' ? window : this);
