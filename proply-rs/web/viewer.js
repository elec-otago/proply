// STEP 3D preview: tessellate the STEP text with occt-import-js (an
// OpenCASCADE WASM build) and render it with three.js.  three.js resolves
// through the import map in index.html; occt-import-js and its .wasm are
// pulled from the CDN on first use.  Any failure leaves the note visible
// with the reason — the download link still works without the preview.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const OCCT_BASE = 'https://cdn.jsdelivr.net/npm/occt-import-js@0.0.23/dist/';

let renderer = null;
let scene = null;
let camera = null;
let controls = null;
let contents = null; // the propeller meshes, replaced per design
let grid = null;
let occt = null;

function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement('script');
    s.src = src;
    s.onload = resolve;
    s.onerror = () => reject(new Error(`cannot load ${src}`));
    document.head.appendChild(s);
  });
}

async function occtInstance() {
  if (occt) {
    return occt;
  }
  if (!window.occtimportjs) {
    await loadScript(`${OCCT_BASE}occt-import-js.js`);
  }
  occt = await window.occtimportjs({
    locateFile: (name) => `${OCCT_BASE}${name}`,
  });
  return occt;
}

function ensureScene(container) {
  if (renderer) {
    return;
  }
  renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(window.devicePixelRatio);
  container.appendChild(renderer.domElement);

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x111318);

  camera = new THREE.PerspectiveCamera(45, 1, 0.1, 1000);
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;

  scene.add(new THREE.HemisphereLight(0xffffff, 0x30343c, 1.4));
  const key = new THREE.DirectionalLight(0xffffff, 2.2);
  key.position.set(1, 2, 1.5);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0x8899bb, 0.8);
  fill.position.set(-1.5, 0.5, -1);
  scene.add(fill);

  contents = new THREE.Group();
  scene.add(contents);

  const resize = () => {
    const w = container.clientWidth;
    const h = container.clientHeight;
    if (w === 0 || h === 0) {
      return;
    }
    renderer.setSize(w, h);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  };
  new ResizeObserver(resize).observe(container);
  resize();
  renderer.setAnimationLoop(() => {
    controls.update();
    renderer.render(scene, camera);
  });
}

/**
 * Remove the on-screen model and grid, leaving the empty scene for the
 * next design.  The renderer stays alive so re-rendering is cheap.
 */
export function clearStep() {
  if (!contents) {
    return;
  }
  for (const child of [...contents.children]) {
    child.geometry.dispose();
    child.material.dispose();
    contents.remove(child);
  }
  if (grid) {
    grid.geometry.dispose();
    scene.remove(grid);
    grid = null;
  }
}

/**
 * Render `stepText` (a STEP document) into `container`; `note` is the
 * overlay element shown while empty or on failure.
 */
export async function showStep(stepText, container, note) {
  try {
    ensureScene(container);

    const importer = await occtInstance();
    const result = importer.ReadStepFile(new TextEncoder().encode(stepText), null);
    if (!result.success || !result.meshes.length) {
      throw new Error('the STEP file contained no readable shapes');
    }

    for (const child of [...contents.children]) {
      child.geometry.dispose();
      child.material.dispose();
      contents.remove(child);
    }

    const box = new THREE.Box3();
    for (const mesh of result.meshes) {
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute(
        'position',
        new THREE.Float32BufferAttribute(mesh.attributes.position.array, 3),
      );
      const hasNormals = !!mesh.attributes.normal;
      if (hasNormals) {
        geometry.setAttribute(
          'normal',
          new THREE.Float32BufferAttribute(mesh.attributes.normal.array, 3),
        );
      }
      geometry.setIndex(new THREE.Uint32BufferAttribute(mesh.index.array, 1));
      const material = new THREE.MeshStandardMaterial({
        color: mesh.color
          ? new THREE.Color(mesh.color[0], mesh.color[1], mesh.color[2])
          : 0xb8bfca,
        metalness: 0.15,
        roughness: 0.55,
        flatShading: !hasNormals,
      });
      const object = new THREE.Mesh(geometry, material);
      contents.add(object);
      box.expandByObject(object);
    }

    // The model is in the scene now — drop the progress overlay.
    note.style.display = 'none';

    // Frame the model: camera on a three-quarter angle, grid underneath.
    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const radius = Math.max(size.length() / 2, 1.0e-6);
    camera.near = radius / 100;
    camera.far = radius * 100;
    // Changing near/far does not rebuild the projection matrix — without
    // this call the camera keeps its initial far=1000, so any prop whose
    // framing distance exceeds 1000 mm-units (roughly radius > 200 mm) is
    // clipped at the far plane and never drawn.
    camera.updateProjectionMatrix();
    camera.position
      .copy(center)
      .add(new THREE.Vector3(1.35, 0.85, 1.7).multiplyScalar(radius * 1.6));
    controls.target.copy(center);
    controls.update();

    if (grid) {
      grid.geometry.dispose();
      scene.remove(grid);
    }
    const step = Math.pow(10, Math.round(Math.log10(radius)));
    grid = new THREE.GridHelper(step * 10, 10, 0x3a4150, 0x232830);
    grid.position.set(center.x, box.min.y - step * 0.02, center.z);
    scene.add(grid);
  } catch (e) {
    note.textContent = `3D preview unavailable: ${e.message || e}`;
    note.style.display = '';
  }
}
