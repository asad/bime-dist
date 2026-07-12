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

    function selectedAtoms(mol) {
        return mol.atoms.filter(function(a) { return !!a.selected; });
    }

    function clearSelection(mol) {
        for (var i = 0; i < mol.atoms.length; i++) mol.atoms[i].selected = false;
        for (var j = 0; j < mol.bonds.length; j++) mol.bonds[j].selected = false;
    }

    function selectedItemCount(mol) {
        var count = 0;
        for (var i = 0; i < mol.atoms.length; i++) if (mol.atoms[i].selected) count++;
        for (var j = 0; j < mol.bonds.length; j++) if (mol.bonds[j].selected) count++;
        return count;
    }

    function translateAtoms(atoms, dx, dy) {
        for (var i = 0; i < atoms.length; i++) {
            atoms[i].x += dx;
            atoms[i].y += dy;
        }
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
        var baseChanged = false;
        if (this._needsNewAtom && !this._startAtom) {
            this.editor.saveHistory();
            var newAtom = mol.addAtom('C', this._needsNewAtom.x, this._needsNewAtom.y);
            this._startAtom = newAtom;
            this._needsNewAtom = null;
            baseChanged = true;
        }
        if (!this._startAtom) return;
        var target = mol.getAtomAt(pos.x, pos.y, 20);
        var endPos;
        if (target && target.id !== this._startAtom.id) {
            endPos = { x: target.x, y: target.y };
        } else {
            endPos = snapAngle(this._startAtom.x, this._startAtom.y, pos.x, pos.y);
        }
        if (baseChanged) this.editor.render();
        this.editor.renderer.drawPreviewBond(this._startAtom.x, this._startAtom.y, endPos.x, endPos.y);
    };

    BondTool.prototype.onMouseUp = function(e, mol, pos) {
        if (!this._dragging) return;
        this._dragging = false;
        this.editor.renderer.clearPreview();
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
        this._dragAtoms = null;
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
            var sel = selectedAtoms(mol);
            if (atom.selected && sel.length > 0) {
                this._dragAtoms = sel;
                this._lastPos = { x: pos.x, y: pos.y };
                this._dragStart = { x: pos.x, y: pos.y };
            } else {
                this._dragAtom = atom;
                this._dragStart = { x: atom.x, y: atom.y };
            }
            this._historySaved = false;
        } else {
            this._panning = true;
            this._lastPos = { x: e.clientX, y: e.clientY };
        }
    };

    MoveTool.prototype.onMouseMove = function(e, mol, pos) {
        if (this._dragAtoms && this._lastPos) {
            var gdx = pos.x - this._lastPos.x;
            var gdy = pos.y - this._lastPos.y;
            if (gdx === 0 && gdy === 0) return;
            if (!this._historySaved && this._dragStart) {
                var tdx = pos.x - this._dragStart.x;
                var tdy = pos.y - this._dragStart.y;
                if (tdx * tdx + tdy * tdy > 0.25) {
                    this.editor.saveHistory();
                    this._historySaved = true;
                }
            }
            translateAtoms(this._dragAtoms, gdx, gdy);
            this._lastPos = { x: pos.x, y: pos.y };
            this.editor.scheduleRender();
        } else if (this._dragAtom) {
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
            this.editor.scheduleRender();
        } else if (this._panning && this._lastPos) {
            var dx = e.clientX - this._lastPos.x;
            var dy = e.clientY - this._lastPos.y;
            this.editor.renderer.offsetX += dx / this.editor.renderer.scale;
            this.editor.renderer.offsetY += dy / this.editor.renderer.scale;
            this._lastPos = { x: e.clientX, y: e.clientY };
            this.editor.scheduleRender();
        }
    };

    MoveTool.prototype.onMouseUp = function(e, mol, pos) {
        // FIX: only fire `changed` if we actually saved history (i.e. the
        // atom genuinely moved). A pure click should not trigger undo state.
        if ((this._dragAtom || this._dragAtoms) && this._historySaved) {
            this.editor.changed();
        }
        this._dragAtom = null;
        this._dragAtoms = null;
        this._panning = false;
        this._lastPos = null;
        this._historySaved = false;
        this._dragStart = null;
    };

    // =========================================================================
    // PanTool — drag the whole molecule (v1.8.9). Discoverable via its own
    // "Pan Mol" button. If atoms are selected, drag that selection as a rigid
    // body; otherwise translate every atom by the cursor delta. This keeps
    // multi-component drawings controllable without losing whole-canvas pan.
    // =========================================================================
    function PanTool(editor) {
        BaseTool.call(this, editor);
        this._dragging = false;
        this._lastPos = null;
        this._dragAtoms = null;
        this._historySaved = false;
        this._dragStart = null;
    }
    PanTool.prototype = Object.create(BaseTool.prototype);

    PanTool.prototype.onMouseDown = function(e, mol, pos) {
        this._dragging = true;
        this._lastPos = { x: pos.x, y: pos.y };
        var sel = selectedAtoms(mol);
        this._dragAtoms = sel.length > 0 ? sel : mol.atoms.slice();
        this._dragStart = { x: pos.x, y: pos.y };
        this._historySaved = false;
    };

    PanTool.prototype.onMouseMove = function(e, mol, pos) {
        if (!this._dragging || !this._lastPos) return;
        var dx = pos.x - this._lastPos.x;
        var dy = pos.y - this._lastPos.y;
        if (dx === 0 && dy === 0) return;
        if (!this._historySaved && this._dragStart) {
            var ddx = pos.x - this._dragStart.x;
            var ddy = pos.y - this._dragStart.y;
            if (ddx * ddx + ddy * ddy > 0.25) {
                this.editor.saveHistory();
                this._historySaved = true;
            }
        }
        translateAtoms(this._dragAtoms || [], dx, dy);
        this._lastPos = { x: pos.x, y: pos.y };
        this.editor.scheduleRender();
    };

    PanTool.prototype.onMouseUp = function(e, mol, pos) {
        if (this._dragging && this._historySaved) {
            this.editor.changed();
        }
        this._dragging = false;
        this._lastPos = null;
        this._dragAtoms = null;
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
            this.editor.scheduleRender();
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
        this.editor.renderer.drawPreviewReactionArrow({
            x1: this._start.x, y1: this._start.y, x2: pos.x, y2: pos.y
        });
    };

    ReactionTool.prototype.onMouseUp = function(e, mol, pos) {
        if (!this._dragging) return;
        this._dragging = false;
        this.editor.renderer.clearPreview();
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
    // Select Tool — click or rectangle selection of atoms and bonds
    // =========================================================================
    function SelectTool(editor) {
        BaseTool.call(this, editor);
        this._startX = 0;
        this._startY = 0;
        this._mode = '';
        this._selectionRect = null;
        this._pressAtom = null;
        this._pressBond = null;
        this._dragAtoms = null;
        this._lastPos = null;
        this._historySaved = false;
        this._additive = false;
        this._touchLike = false;
    }
    SelectTool.prototype = Object.create(BaseTool.prototype);

    SelectTool.prototype.deactivate = function() {
        if (this._selectionRect) {
            this._selectionRect.remove();
            this._selectionRect = null;
        }
        this._mode = '';
        this._pressAtom = null;
        this._pressBond = null;
        this._dragAtoms = null;
        this._lastPos = null;
        this._historySaved = false;
        this.active = false;
    };

    SelectTool.prototype._hitRadii = function(e) {
        var scale = this.editor && this.editor.renderer ? (this.editor.renderer.scale || 1) : 1;
        var touch = !!(e && ((e.touches && e.touches.length) ||
            (e.changedTouches && e.changedTouches.length) || e.pointerType === 'touch'));
        return {
            atom: touch ? Math.max(10, 18 / scale) : Math.max(6, 12 / scale),
            bond: touch ? Math.max(7, 12 / scale) : Math.max(4, 7 / scale),
            drag: touch ? Math.max(4, 7 / scale) : Math.max(2, 4 / scale)
        };
    };

    SelectTool.prototype._createSelectionRect = function() {
        if (this._selectionRect) return;
        var svg = this.editor.renderer.svg;
        this._selectionRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        this._selectionRect.setAttribute('fill', 'rgba(13,148,136,0.1)');
        this._selectionRect.setAttribute('stroke', '#0d9488');
        this._selectionRect.setAttribute('stroke-width', '1');
        this._selectionRect.setAttribute('stroke-dasharray', '4,2');
        svg.appendChild(this._selectionRect);
    };

    SelectTool.prototype._updateSelectionRect = function(pos) {
        if (!this._selectionRect) return;
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

    SelectTool.prototype.onMouseDown = function(e, mol, pos) {
        this._startX = pos.x;
        this._startY = pos.y;
        this._mode = 'press';
        this._historySaved = false;
        this._lastPos = { x: pos.x, y: pos.y };
        this._additive = !!(e && (e.shiftKey || e.ctrlKey || e.metaKey));
        this._touchLike = !!(e && ((e.touches && e.touches.length) ||
            (e.changedTouches && e.changedTouches.length) || e.pointerType === 'touch'));
        var radii = this._hitRadii(e);
        this._pressAtom = mol.getAtomAt(pos.x, pos.y, radii.atom);
        this._pressBond = this._pressAtom ? null : mol.getBondAt(pos.x, pos.y, radii.bond);
        this._dragAtoms = null;

        if (this._pressAtom) {
            if (!this._additive && !this._pressAtom.selected) {
                clearSelection(mol);
                this._pressAtom.selected = true;
            }
            var sel = selectedAtoms(mol);
            this._dragAtoms = (this._pressAtom.selected && sel.length > 0) ? sel : [this._pressAtom];
        } else if (this._pressBond) {
            if (!this._additive && !this._pressBond.selected) {
                clearSelection(mol);
                this._pressBond.selected = true;
            }
            var a1 = mol.getAtom(this._pressBond.atom1);
            var a2 = mol.getAtom(this._pressBond.atom2);
            var selected = selectedAtoms(mol);
            this._dragAtoms = (this._pressBond.selected && selected.length > 0)
                ? selected
                : (a1 && a2 ? [a1, a2] : []);
        }
    };

    SelectTool.prototype.onMouseMove = function(e, mol, pos) {
        if (!this._mode) return;
        var dx0 = pos.x - this._startX;
        var dy0 = pos.y - this._startY;
        var radii = this._hitRadii(e);

        if (this._mode === 'press') {
            if (dx0 * dx0 + dy0 * dy0 <= radii.drag * radii.drag) return;
            if (this._dragAtoms && this._dragAtoms.length > 0) {
                this._mode = 'drag';
                this.editor.saveHistory();
                this._historySaved = true;
            } else {
                this._mode = 'box';
                this._createSelectionRect();
            }
        }

        if (this._mode === 'drag' && this._lastPos && this._dragAtoms) {
            var dx = pos.x - this._lastPos.x;
            var dy = pos.y - this._lastPos.y;
            if (dx === 0 && dy === 0) return;
            translateAtoms(this._dragAtoms, dx, dy);
            this._lastPos = { x: pos.x, y: pos.y };
            this.editor.scheduleRender();
        } else if (this._mode === 'box') {
            this._updateSelectionRect(pos);
        }
    };

    SelectTool.prototype.onMouseUp = function(e, mol, pos) {
        if (!this._mode) return;

        var dx = pos.x - this._startX;
        var dy = pos.y - this._startY;
        var scale = this.editor && this.editor.renderer ? (this.editor.renderer.scale || 1) : 1;
        var clickTol = this._touchLike ? Math.max(4, 7 / scale) : Math.max(3, 5 / scale);

        if (this._mode === 'drag') {
            if (this._historySaved) {
                this.editor.changed();
            }
            this.editor.showInfo(this._dragAtoms.length + ' atom' +
                (this._dragAtoms.length === 1 ? '' : 's') + ' moved');
            this._mode = '';
            this._pressAtom = null;
            this._pressBond = null;
            this._dragAtoms = null;
            this._lastPos = null;
            this._historySaved = false;
            return;
        }

        if (this._mode === 'press' && dx * dx + dy * dy <= clickTol * clickTol) {
            if (!this._additive && !this._pressAtom && !this._pressBond) {
                clearSelection(mol);
            } else if (this._additive) {
                if (this._pressAtom) {
                    this._pressAtom.selected = !this._pressAtom.selected;
                } else if (this._pressBond) {
                    this._pressBond.selected = !this._pressBond.selected;
                }
            } else if (this._pressAtom) {
                clearSelection(mol);
                this._pressAtom.selected = true;
            } else if (this._pressBond) {
                clearSelection(mol);
                this._pressBond.selected = true;
            }
            if (this._selectionRect) {
                this._selectionRect.remove();
                this._selectionRect = null;
            }
            var itemCount = selectedItemCount(mol);
            this.editor.render();
            this.editor.showInfo(itemCount + ' item' + (itemCount === 1 ? '' : 's') + ' selected');
            this._mode = '';
            this._pressAtom = null;
            this._pressBond = null;
            this._dragAtoms = null;
            this._lastPos = null;
            return;
        }

        // Find atoms within the selection rectangle (in molecule coords)
        var minX = Math.min(this._startX, pos.x);
        var maxX = Math.max(this._startX, pos.x);
        var minY = Math.min(this._startY, pos.y);
        var maxY = Math.max(this._startY, pos.y);

        // Select atoms in range
        var selectedCount = 0;
        if (!this._additive) clearSelection(mol);
        for (var i = 0; i < mol.atoms.length; i++) {
            var a = mol.atoms[i];
            if (a.x >= minX && a.x <= maxX && a.y >= minY && a.y <= maxY) {
                a.selected = true;
                selectedCount++;
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
        this._mode = '';
        this._pressAtom = null;
        this._pressBond = null;
        this._dragAtoms = null;
        this._lastPos = null;
    };

    // =========================================================================
    // v2.0.29: CurlyArrowTool — IUPAC mechanism curly arrows.
    // Click first atom to set the source, click second atom to complete the
    // arrow. Press Escape (or click again on the source) to cancel. Style
    // ('pair' or 'single') is passed at construction so the same tool class
    // handles both with no branching in the editor's dispatch.
    // =========================================================================

    function CurlyArrowTool(editor, style) {
        BaseTool.call(this, editor);
        this.style = (style === 'single') ? 'single' : 'pair';
        this._pendingSource = null;     // endpoint ref (atom or bond) that started the arrow
        this._pendingHighlight = null;  // the live object (atom or bond) carrying the flag
    }
    CurlyArrowTool.prototype = Object.create(BaseTool.prototype);

    CurlyArrowTool.prototype.activate = function() {
        this._pendingSource = null;
        this._pendingHighlight = null;
    };

    CurlyArrowTool.prototype.deactivate = function() {
        this._clearHighlight();
        this._pendingSource = null;
    };

    // v2.0.31: hit-test atoms first, then bonds; both can be arrow endpoints.
    CurlyArrowTool.prototype._pickEndpoint = function(mol, pos) {
        var atom = mol.getAtomAt(pos.x, pos.y);
        if (atom) {
            return { ref: { kind: 'atom', atomId: atom.id }, live: atom };
        }
        if (typeof mol.getBondAt === 'function') {
            var bond = mol.getBondAt(pos.x, pos.y);
            if (bond) {
                return { ref: { kind: 'bond', bondId: bond.id }, live: bond };
            }
        }
        return null;
    };

    function sameRef(a, b) {
        if (!a || !b) { return false; }
        if (a.kind !== b.kind) { return false; }
        if (a.kind === 'atom') { return a.atomId === b.atomId; }
        return a.bondId === b.bondId;
    }

    CurlyArrowTool.prototype.onMouseDown = function(e, mol, pos) {
        var picked = this._pickEndpoint(mol, pos);
        if (!picked) {
            // Click on empty canvas cancels the in-flight arrow.
            this._clearHighlight();
            this._pendingSource = null;
            this.editor.render();
            return;
        }
        if (!this._pendingSource) {
            this._pendingSource = picked.ref;
            this._pendingHighlight = picked.live;
            picked.live._curlyArrowPending = true;
            this.editor.render();
            return;
        }
        if (sameRef(picked.ref, this._pendingSource)) {
            // Re-click on source = cancel.
            this._clearHighlight();
            this._pendingSource = null;
            this.editor.render();
            return;
        }
        // Two distinct endpoints picked → save history, add arrow, reset.
        this.editor.saveHistory();
        var fromRef = this._pendingSource;
        this._clearHighlight();
        this._pendingSource = null;
        try {
            mol.addCurlyArrow({
                style: this.style,
                from: fromRef,
                to:   picked.ref
            });
        } catch (err) {
            if (typeof console !== 'undefined' && console.warn) {
                console.warn('CurlyArrowTool: addCurlyArrow rejected:', err && err.message);
            }
        }
        this.editor.changed();
    };

    CurlyArrowTool.prototype._clearHighlight = function() {
        if (this._pendingHighlight) { this._pendingHighlight._curlyArrowPending = false; }
        this._pendingHighlight = null;
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
        PanTool: PanTool,
        RingTool: RingTool,
        ChainTool: ChainTool,
        StereoBondTool: StereoBondTool,
        ReactionTool: ReactionTool,
        MapTool: MapTool,
        SelectTool: SelectTool,
        CurlyArrowTool: CurlyArrowTool
    };

})(typeof window !== 'undefined' ? window : this);
