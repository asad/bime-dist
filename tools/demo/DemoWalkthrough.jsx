const demoPoster = 'screenshots/bime-demo-poster.png';

export const demoScenes = [
  {
    id: 'start-fast',
    image: demoPoster,
    sourceFrame: '01-workbench.png',
    eyebrow: '01',
    title: 'Start Fast',
    body: 'Draw, edit, clean, zoom, and inspect molecules in one focused workbench.',
    highlights: ['Atoms', 'Bonds', 'Clean 2D', 'Inspector'],
  },
  {
    id: 'properties',
    image: demoPoster,
    sourceFrame: '02-properties.png',
    eyebrow: '02',
    title: 'Properties',
    body: 'Name structures, add export comments, and decide whether labels appear on canvas.',
    highlights: ['Name', 'Comment', 'Show label', 'MOL/RXN metadata'],
  },
  {
    id: 'quick-load',
    image: demoPoster,
    sourceFrame: '03-quick-load.png',
    eyebrow: '03',
    title: 'Quick Load',
    body: 'Open common molecules and reaction examples without leaving the editor.',
    highlights: ['Molecules', 'Reactions', 'One-click load'],
  },
  {
    id: 'library-search',
    image: demoPoster,
    sourceFrame: '04-library-search.png',
    eyebrow: '04',
    title: 'Library Search',
    body: 'Search the built-in molecule library and load useful hits quickly.',
    highlights: ['Exact', 'Substructure', 'Similarity', 'MCS'],
  },
  {
    id: 'reactions',
    image: demoPoster,
    sourceFrame: '05-reaction-properties.png',
    eyebrow: '05',
    title: 'Reactions',
    body: 'Add reaction conditions, yield, and step notes as editable properties.',
    highlights: ['Conditions', 'Yield', 'Step note', 'Arrow label'],
  },
  {
    id: 'export',
    image: demoPoster,
    sourceFrame: '06-image-export.png',
    eyebrow: '06',
    title: 'Export',
    body: 'Preview images and switch between SMILES, MOL/RXN, SDF, and image outputs.',
    highlights: ['SMILES', 'MOL/RXN', 'SDF', 'Image'],
  },
  {
    id: 'canvas-mode',
    image: demoPoster,
    sourceFrame: '07-pathway-top.png',
    eyebrow: '07',
    title: 'Canvas Mode',
    body: 'Toggle into reaction and pathway mode when the work grows beyond one molecule.',
    highlights: ['Pathway toggle', 'Reaction map', 'Extended canvas'],
  },
  {
    id: 'pathway-drafting',
    image: demoPoster,
    sourceFrame: '09-pathway-canvas-deep.png',
    eyebrow: '08',
    title: 'Pathway Drafting',
    body: 'Use structures, residues, arrows, steps, notes, and layout controls together.',
    highlights: ['Residues', 'Arrows', 'Steps', 'Notes', 'Auto layout'],
  },
];

export const demoVoiceover = [
  'BIME opens directly into the Workbench, so you can start drawing immediately.',
  'Use the focused toolbar for atoms, bonds, rings, edit actions, clean layout, zoom, and fit.',
  'The Properties tab keeps identity metadata in one place.',
  'Name the molecule or reaction, add an export comment, and choose whether the name appears on the canvas.',
  'Quick Load gives you ready molecule and reaction examples.',
  'The built-in library search helps you find structures and load them into the editor without leaving the page.',
  'For reactions, add conditions, yield, and step notes as editable properties.',
  'Preview images and export SMILES, MOL/RXN, SDF, SVG, and PNG.',
  'When you need a bigger drawing surface, switch into Reaction and Pathway Canvas mode.',
  'Build maps with structures, residues, arrows, mechanism steps, compartments, notes, import references, auto layout, and SVG export.',
  'BIME is designed to stay fast, local, and clear across desktop, tablet, and mobile.',
];

const defaultBrand = {
  product: 'BIME Workbench',
  tagline: 'A fast local molecule, reaction, and pathway editor',
  cta: 'Try the Workbench',
};

export function DemoWalkthrough({
  scenes = demoScenes,
  brand = defaultBrand,
  activeScene = 0,
  onSceneSelect,
  showTimeline = true,
  showVoiceover = false,
}) {
  const safeScenes = Array.isArray(scenes) && scenes.length ? scenes : demoScenes;
  const safeSceneIndex = Math.max(0, Math.min(activeScene, safeScenes.length - 1));
  const scene = safeScenes[safeSceneIndex];

  return (
    <article className="demo-walkthrough" aria-label={`${brand.product} demo walkthrough`}>
      <header className="demo-hero">
        <p className="demo-kicker">User demo</p>
        <h1>{brand.product}</h1>
        <p>{brand.tagline}</p>
      </header>

      <section className="demo-stage" aria-labelledby={`demo-scene-${scene.id}`}>
        <img className="demo-stage-image" src={scene.image} alt={`${scene.title} screen`} />
        <div className="demo-caption">
          <span className="demo-number">{scene.eyebrow}</span>
          <div>
            <h2 id={`demo-scene-${scene.id}`}>{scene.title}</h2>
            <p>{scene.body}</p>
            <ul className="demo-highlight-list" aria-label={`${scene.title} highlights`}>
              {scene.highlights.map((highlight) => (
                <li key={highlight}>{highlight}</li>
              ))}
            </ul>
          </div>
        </div>
      </section>

      {showTimeline ? (
        <nav className="demo-timeline" aria-label="Demo scenes">
          {safeScenes.map((item, index) => (
            <button
              type="button"
              key={item.id}
              className={index === safeSceneIndex ? 'is-active' : undefined}
              aria-current={index === safeSceneIndex ? 'step' : undefined}
              onClick={() => {
                if (onSceneSelect) onSceneSelect(index);
              }}
            >
              <span>{item.eyebrow}</span>
              {item.title}
            </button>
          ))}
        </nav>
      ) : null}

      {showVoiceover ? (
        <section className="demo-voiceover" aria-label="Voiceover script">
          {demoVoiceover.map((line) => (
            <p key={line}>{line}</p>
          ))}
        </section>
      ) : null}

      <footer className="demo-footer">
        <a href="workbench.html">{brand.cta}</a>
      </footer>
    </article>
  );
}

export default DemoWalkthrough;
