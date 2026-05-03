/**
 * Tools.js — Interactive drawing tools for molecule editor
 * * Copyright (c) 2026 BioInception PVT LTD, Cambridge, UK and Syed Asad Rahman.
 * All rights reserved. Licensed under the Apache License, Version 2.0 — see LICENSE.txt
 * If you use BIME, please cite: S. A. Rahman, BIME: BioInception Molecular Editor (2026)
 *
 * Each tool is a state machine handling mousedown/mousemove/mouseup events.
 */
(function(global) {
    'use strict';

    var BOND_LENGTH = Molecule.BOND_LENGTH;
    var SNAP_ANGLES = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330];

    // Snap a point to the nearest 30-degree angle from an origin
    function snapAngle(ox, oy, x, y) {
        var dx = x - ox, dy = y - oy;
        var dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < 5) return { x: x, y: y };
        var angle = Math.atan2(dy, dx) * 180 / Math.PI;
        var best = SNAP_ANGLES[0], bestDiff = 999;
        for (var i = 0; i < SNAP_ANGLES.length; i++) {
            var diff = Math.abs(((angle - SNAP_ANGLES[i] + 540) % 360) - 180);
            if (diff < bestDiff) { bestDiff = diff; best = SNAP_ANGLES[i]; }
        }
        var rad = best * Math.PI / 180;
        return { x: ox + BOND_LENGTH * Math.cos(rad), y: oy + BOND_LENGTH * Math.sin(rad) };
    }

    // =========================================================================
    // Base Tool
    // =========================================================================
    function BaseTool(editor) {
        this.editor = editor;
        this.active = false;
    }
    BaseTool.prototype.onMouseDown = function(e, mol, pos) {};
    BaseTool.prototype.onMouseMove = function(e, mol, pos) {};
    BaseTool.prototype.onMouseUp = function(e, mol, pos) {};
    BaseTool.prototype.onKeyDown = function(e, mol) {};
    BaseTool.prototype.activate = function() { this.active = true; };
    BaseTool.prototype.deactivate = function() { this.active = false; };

    // =========================================================================
    // Atom Tool — place atoms, change element of existing atoms
    // =========================================================================
    function AtomTool(editor, symbol) {
        BaseTool.call(this, editor);
        this.symbol = symbol || 'C';
    }
    AtomTool.prototype = Object.create(BaseTool.prototype);

    AtomTool.prototype.onMouseDown = function(e, mol, pos) {
        var existing = mol.getAtomAt(pos.x, pos.y);
        if (existing) {
            // Change element of existing atom
            if (existing.symbol !== this.symbol) {
                this.editor.saveHistory();
                existing.symbol = this.symbol;
                this.editor.changed();
            }
        } else {
            // Place new atom
            this.editor.saveHistory();
            mol.addAtom(this.symbol, pos.x, pos.y);
            this.editor.changed();
        }
    };

    // =========================================================================
    // Bond Tool — draw bonds between atoms
    // =========================================================================
    function BondTool(editor, bondType) {
        BaseTool.call(this, editor);
        this.bondType = bondType || Molecule.BOND_SINGLE;
        this._startAtom = null;
        this._dragging = false;
    }
    BondTool.prototype = Object.create(BaseTool.prototype);

    BondTool.prototype.onMouseDown = function(e, mol, pos) {
        var atom = mol.getAtomAt(pos.x, pos.y);
        if (atom) {
            this._startAtom = atom;
            this._dragging = true;
        } else {
            // FIX: don't save history here; onMouseUp will save it once for the entire operation
            // (previously this caused double undo entries — one for atom creation, one for bond)
            this._needsNewAtom = { x: pos.x, y: pos.y };
            this._startAtom = null;
            this._dragging = true;
        }
    };

    BondTool.prototype.onMouseMove = function(e, mol, pos) {
        if (!this._dragging) return;
        // FIX: create deferred start atom on first drag move if needed
        if (this._needsNewAtom && !this._startAtom) {
            this.editor.saveHistory();
            var newAtom = mol.addAtom('C', this._needsNewAtom.x, this._needsNewAtom.y);
            this._startAtom = newAtom;
            this._needsNewAtom = null;
        }
        if (!this._startAtom) return;
        var target = mol.getAtomAt(pos.x, pos.y, 20);
        var endPos;
        if (target && target.id !== this._startAtom.id) {
            endPos = { x: target.x, y: target.y };
        } else {
            endPos = snapAngle(this._startAtom.x, this._startAtom.y, pos.x, pos.y);
        }
        this.editor.renderer.render();
        this.editor.renderer.drawPreviewBond(this._startAtom.x, this._startAtom.y, endPos.x, endPos.y);
    };

    BondTool.prototype.onMouseUp = function(e, mol, pos) {
        if (!this._dragging) return;
        this._dragging = false;
        // FIX: create deferred start atom if mouse was released without moving
        if (this._needsNewAtom && !this._startAtom) {
            this.editor.saveHistory();
            var newAtom = mol.addAtom('C', this._needsNewAtom.x, this._needsNewAtom.y);
            this._startAtom = newAtom;
            this._needsNewAtom = null;
            this.editor.changed();
            this._startAtom = null;
            return;
        }
        if (!this._startAtom) return;

        var target = mol.getAtomAt(pos.x, pos.y, 20);
        if (target && target.id === this._startAtom.id) {
            // Clicked same atom — try to cycle bond type on an existing bond
            var bonds = mol.getBondsOfAtom(this._startAtom.id);
            if (bonds.length === 1) {
                this.editor.saveHistory();
                bonds[0].type = bonds[0].type % 3 + 1;
                this.editor.changed();
            }
            this._startAtom = null;
            return;
        }

        this.editor.saveHistory();

        var endAtom;
        if (target) {
            endAtom = target;
        } else {
            var endPos = snapAngle(this._startAtom.x, this._startAtom.y, pos.x, pos.y);
            endAtom = mol.addAtom('C', endPos.x, endPos.y);
        }

        var existingBond = mol.getBondBetween(this._startAtom.id, endAtom.id);
        if (existingBond) {
            // Cycle bond type
            existingBond.type = existingBond.type % 3 + 1;
        } else {
            mol.addBond(this._startAtom.id, endAtom.id, this.bondType);
        }

        this._startAtom = null;
        this.editor.changed();
    };

    // =========================================================================
    // Delete Tool — click atom or bond to delete
    // =========================================================================
    function DeleteTool(editor) {
        BaseTool.call(this, editor);
    }
    DeleteTool.prototype = Object.create(BaseTool.prototype);

    DeleteTool.prototype.onMouseDown = function(e, mol, pos) {
        var atom = mol.getAtomAt(pos.x, pos.y);
        if (atom) {
            this.editor.saveHistory();
            mol.removeAtom(atom.id);
            this.editor.changed();
            return;
        }
        var bond = mol.getBondAt(pos.x, pos.y);
        if (bond) {
            this.editor.saveHistory();
            mol.removeBond(bond.id);
            this.editor.changed();
        }
    };

    // =========================================================================
    // Move Tool — drag atoms, pan canvas
    // =========================================================================
    function MoveTool(editor) {
        BaseTool.call(this, editor);
        this._dragAtom = null;
        this._panning = false;
        this._lastPos = null;
        this._historySaved = false;
        this._dragStart = null;
    }
    MoveTool.prototype = Object.create(BaseTool.prototype);

    MoveTool.prototype.onMouseDown = function(e, mol, pos) {
        var atom = mol.getAtomAt(pos.x, pos.y);
        if (atom) {
            // FIX: defer saveHistory until the user actually moves the atom.
            // Previously a click-and-release with no movement still pushed an
            // identical snapshot to history, causing a no-op undo entry.
            this._dragAtom = atom;
            this._dragStart = { x: atom.x, y: atom.y };
            this._historySaved = false;
        } else {
            this._panning = true;
            this._lastPos = { x: e.clientX, y: e.clientY };
        }
    };

    MoveTool.prototype.onMouseMove = function(e, mol, pos) {
        if (this._dragAtom) {
            // FIX: only save history once, and only when the atom actually moves.
            // Use small tolerance instead of strict float equality.
            if (!this._historySaved && this._dragStart) {
                var ddx = pos.x - this._dragStart.x;
                var ddy = pos.y - this._dragStart.y;
                if (ddx * ddx + ddy * ddy > 0.25) {
                    this.editor.saveHistory();
                    this._historySaved = true;
                }
            }
            this._dragAtom.x = pos.x;
            this._dragAtom.y = pos.y;
            this.editor.render();
        } else if (this._panning && this._lastPos) {
            var dx = e.clientX - this._lastPos.x;
            var dy = e.clientY - this._lastPos.y;
            this.editor.renderer.offsetX += dx / this.editor.renderer.scale;
            this.editor.renderer.offsetY += dy / this.editor.renderer.scale;
            this._lastPos = { x: e.clientX, y: e.clientY };
            this.editor.render();
        }
    };

    MoveTool.prototype.onMouseUp = function(e, mol, pos) {
        // FIX: only fire `changed` if we actually saved history (i.e. the
        // atom genuinely moved). A pure click should not trigger undo state.
        if (this._dragAtom && this._historySaved) {
            this.editor.changed();
        }
        this._dragAtom = null;
        this._panning = false;
        this._lastPos = null;
        this._historySaved = false;
        this._dragStart = null;
    };

    // =========================================================================
    // Ring Tool — place ring templates
    // =========================================================================
    function RingTool(editor, size) {
        BaseTool.call(this, editor);
        this.ringSize = size || 6;
    }
    RingTool.prototype = Object.create(BaseTool.prototype);

    RingTool.prototype.onMouseDown = function(e, mol, pos) {
        var atom = mol.getAtomAt(pos.x, pos.y);
        this.editor.saveHistory();
        if (atom) {
            mol.addRing(this.ringSize, 0, 0, atom.id);
        } else {
            mol.addRing(this.ringSize, pos.x, pos.y);
        }
        this.editor.changed();
    };

    // =========================================================================
    // Chain Tool — drag to draw carbon chains
    // =========================================================================
    function ChainTool(editor) {
        BaseTool.call(this, editor);
        this._atoms = [];
        this._dragging = false;
    }
    ChainTool.prototype = Object.create(BaseTool.prototype);

    ChainTool.prototype.onMouseDown = function(e, mol, pos) {
        // FIX: defer saveHistory until we actually modify the molecule (onMouseMove or onMouseUp)
        this._historySaved = false;
        var atom = mol.getAtomAt(pos.x, pos.y);
        if (!atom) {
            this.editor.saveHistory();
            this._historySaved = true;
            atom = mol.addAtom('C', pos.x, pos.y);
        }
        this._atoms = [atom];
        this._dragging = true;
    };

    ChainTool.prototype.onMouseMove = function(e, mol, pos) {
        if (!this._dragging || this._atoms.length === 0) return;
        var last = this._atoms[this._atoms.length - 1];
        var dx = pos.x - last.x, dy = pos.y - last.y;
        var dist = Math.sqrt(dx * dx + dy * dy);
        if (dist >= BOND_LENGTH * 0.8) {
            // Save history on first actual chain extension
            if (!this._historySaved) {
                this.editor.saveHistory();
                this._historySaved = true;
            }
            var snapped = snapAngle(last.x, last.y, pos.x, pos.y);
            var existing = mol.getAtomAt(snapped.x, snapped.y, 10);
            var next = existing || mol.addAtom('C', snapped.x, snapped.y);
            if (!mol.getBondBetween(last.id, next.id)) {
                mol.addBond(last.id, next.id, Molecule.BOND_SINGLE);
            }
            this._atoms.push(next);
            this.editor.render();
        }
    };

    ChainTool.prototype.onMouseUp = function(e, mol, pos) {
        this._dragging = false;
        this._atoms = [];
        this.editor.changed();
    };

    // =========================================================================
    // Stereo Bond Tool — cycle stereo on existing bonds
    // =========================================================================
    function StereoBondTool(editor) {
        BaseTool.call(this, editor);
    }
    StereoBondTool.prototype = Object.create(BaseTool.prototype);

    StereoBondTool.prototype.onMouseDown = function(e, mol, pos) {
        var bond = mol.getBondAt(pos.x, pos.y);
        if (bond) {
            this.editor.saveHistory();
            if (bond.stereo === Molecule.STEREO_NONE) bond.stereo = Molecule.STEREO_WEDGE;
            else if (bond.stereo === Molecule.STEREO_WEDGE) bond.stereo = Molecule.STEREO_DASH;
            else bond.stereo = Molecule.STEREO_NONE;
            this.editor.changed();
        }
    };

    // =========================================================================
    // Reaction Arrow Tool — draw reaction arrow
    // =========================================================================
    function ReactionTool(editor) {
        BaseTool.call(this, editor);
        this._start = null;
        this._dragging = false;
    }
    ReactionTool.prototype = Object.create(BaseTool.prototype);

    ReactionTool.prototype.onMouseDown = function(e, mol, pos) {
        this._start = { x: pos.x, y: pos.y };
        this._dragging = true;
    };

    ReactionTool.prototype.onMouseMove = function(e, mol, pos) {
        if (!this._dragging) return;
        this.editor.renderer.render();
        this.editor.renderer._drawReactionArrow({
            x1: this._start.x, y1: this._start.y, x2: pos.x, y2: pos.y
        });
    };

    ReactionTool.prototype.onMouseUp = function(e, mol, pos) {
        if (!this._dragging) return;
        this._dragging = false;
        var dx = pos.x - this._start.x;
        if (Math.abs(dx) > 20) {
            this.editor.saveHistory();
            mol.reactionArrow = { x1: this._start.x, y1: this._start.y, x2: pos.x, y2: pos.y };
            this.editor.changed();
        }
        this._start = null;
    };

    // =========================================================================
    // Atom-Atom Mapping Tool
    // =========================================================================
    function MapTool(editor) {
        BaseTool.call(this, editor);
        this._nextMap = 1;
        this._firstAtom = null;
    }
    MapTool.prototype = Object.create(BaseTool.prototype);

    MapTool.prototype.onMouseDown = function(e, mol, pos) {
        var atom = mol.getAtomAt(pos.x, pos.y);
        if (!atom) return;

        if (!this._firstAtom) {
            // First click — select reactant atom
            this._firstAtom = atom;
            atom.highlighted = true;
            this.editor.render();
        } else {
            // Second click — map the pair
            this.editor.saveHistory();
            this._firstAtom.mapNumber = this._nextMap;
            atom.mapNumber = this._nextMap;
            this._nextMap++;
            this._firstAtom.highlighted = false;
            this._firstAtom = null;
            this.editor.changed();
        }
    };

    MapTool.prototype.deactivate = function() {
        if (this._firstAtom) {
            this._firstAtom.highlighted = false;
            this._firstAtom = null;
        }
        this.active = false;
    };

    // =========================================================================
    // Select Tool — rectangle selection of atoms and bonds
    // =========================================================================
    function SelectTool(editor) {
        BaseTool.call(this, editor);
        this._startX = 0;
        this._startY = 0;
        this._selecting = false;
        this._selectionRect = null;
    }
    SelectTool.prototype = Object.create(BaseTool.prototype);

    SelectTool.prototype.deactivate = function() {
        if (this._selectionRect) {
            this._selectionRect.remove();
            this._selectionRect = null;
        }
        this.active = false;
    };

    SelectTool.prototype.onMouseDown = function(e, mol, pos) {
        this._startX = pos.x;
        this._startY = pos.y;
        this._selecting = true;
        // Create visual selection rectangle
        var svg = this.editor.renderer.svg;
        this._selectionRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        this._selectionRect.setAttribute('fill', 'rgba(13,148,136,0.1)');
        this._selectionRect.setAttribute('stroke', '#0d9488');
        this._selectionRect.setAttribute('stroke-width', '1');
        this._selectionRect.setAttribute('stroke-dasharray', '4,2');
        svg.appendChild(this._selectionRect);
    };

    SelectTool.prototype.onMouseMove = function(e, mol, pos) {
        if (!this._selecting || !this._selectionRect) return;
        var renderer = this.editor.renderer;
        // Convert mol coords to screen coords for the rectangle overlay
        var sx1 = (this._startX + renderer.offsetX) * renderer.scale;
        var sy1 = (this._startY + renderer.offsetY) * renderer.scale;
        var sx2 = (pos.x + renderer.offsetX) * renderer.scale;
        var sy2 = (pos.y + renderer.offsetY) * renderer.scale;
        var minX = Math.min(sx1, sx2);
        var minY = Math.min(sy1, sy2);
        var w = Math.abs(sx2 - sx1);
        var h = Math.abs(sy2 - sy1);
        this._selectionRect.setAttribute('x', minX);
        this._selectionRect.setAttribute('y', minY);
        this._selectionRect.setAttribute('width', w);
        this._selectionRect.setAttribute('height', h);
    };

    SelectTool.prototype.onMouseUp = function(e, mol, pos) {
        if (!this._selecting) return;
        this._selecting = false;

        // Find atoms within the selection rectangle (in molecule coords)
        var minX = Math.min(this._startX, pos.x);
        var maxX = Math.max(this._startX, pos.x);
        var minY = Math.min(this._startY, pos.y);
        var maxY = Math.max(this._startY, pos.y);

        // Select atoms in range
        var selectedCount = 0;
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            if (a.x >= minX && a.x <= maxX && a.y >= minY && a.y <= maxY) {
                a.selected = true;
                selectedCount++;
            } else {
                a.selected = false;
            }
        }
        // Select bonds where both atoms are selected
        for (var j = 0; j < mol.bonds.length; j++) {
            var b = mol.bonds[j];
            var a1 = mol.getAtom(b.atom1);
            var a2 = mol.getAtom(b.atom2);
            b.selected = !!(a1 && a1.selected && a2 && a2.selected);
        }

        // Remove selection rect
        if (this._selectionRect) {
            this._selectionRect.remove();
            this._selectionRect = null;
        }

        this.editor.render();
        this.editor.showInfo(selectedCount + ' atoms selected');
    };

    // =========================================================================
    // Exports
    // =========================================================================
    global.EditorTools = {
        BaseTool: BaseTool,
        AtomTool: AtomTool,
        BondTool: BondTool,
        DeleteTool: DeleteTool,
        MoveTool: MoveTool,
        RingTool: RingTool,
        ChainTool: ChainTool,
        StereoBondTool: StereoBondTool,
        ReactionTool: ReactionTool,
        MapTool: MapTool,
        SelectTool: SelectTool
    };

})(typeof window !== 'undefined' ? window : this);
