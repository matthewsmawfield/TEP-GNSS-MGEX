#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

async function buildStaticSite() {
    console.log('🔨 Building static site...');

    try {
        const distDir = path.join(__dirname, 'dist');
        if (fs.existsSync(distDir)) {
            console.log('🧹 Cleaning dist directory...');
            fs.rmSync(distDir, { recursive: true, force: true });
        }

        const manifestPath = path.join(__dirname, 'manifest.json');
        const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));

        const indexPath = path.join(__dirname, 'index.html');
        let indexContent = fs.readFileSync(indexPath, 'utf8');

        let componentsHtml = '';
        const sortedSections = manifest.sections.sort((a, b) => a.order - b.order);

        for (const section of sortedSections) {
            console.log(`📄 Loading section: ${section.title}`);
            const componentPath = path.join(__dirname, 'components', section.file);
            if (fs.existsSync(componentPath)) {
                const componentHtml = fs.readFileSync(componentPath, 'utf8');
                componentsHtml += `
                <section id="${section.id}" class="manuscript-section" data-section="${section.title}">
                    ${componentHtml}
                </section>`;
            } else {
                console.warn(`⚠️  Component not found: ${section.file}`);
                componentsHtml += `
                <section id="${section.id}" class="manuscript-section" data-section="${section.title}">
                    <div style="background-color: #ffe6e6; border: 1px solid #ff9999; padding: 15px; margin: 20px 0; border-radius: 5px;">
                        <h3 style="color: #cc0000; margin-top: 0;">Missing Section: ${section.title}</h3>
                        <p>Component file <code>components/${section.file}</code> not found</p>
                    </div>
                </section>`;
            }
        }

        const staticContent = indexContent
            .replace(
                /<div id="loading".*?<\/div>\s*<div id="manuscript-content".*?<\/div>/s,
                `<div id="manuscript-content">${componentsHtml}</div>`
            )
            .replace(
                /<!-- Component Loading Script -->[\s\S]*?<\/script>/,
                '<!-- Static build - components pre-loaded -->'
            )
            .replace(
                '<main id="main-content" role="main">',
                '<!-- This is a statically built version for SEO/deployment -->\n    <main id="main-content" role="main">'
            );

        if (!fs.existsSync(distDir)) {
            fs.mkdirSync(distDir, { recursive: true });
        }

        const outputPath = path.join(distDir, 'index.html');
        fs.writeFileSync(outputPath, staticContent, 'utf8');
        console.log(`✅ Built: ${outputPath}`);

        // Copy public assets
        const publicDir = path.join(__dirname, 'public');
        if (fs.existsSync(publicDir)) {
            const distPublicDir = path.join(distDir, 'public');
            const copyRecursive = (src, dest) => {
                if (!fs.existsSync(dest)) fs.mkdirSync(dest, { recursive: true });
                for (const entry of fs.readdirSync(src)) {
                    const srcPath = path.join(src, entry);
                    const destPath = path.join(dest, entry);
                    if (fs.statSync(srcPath).isDirectory()) {
                        copyRecursive(srcPath, destPath);
                    } else {
                        fs.copyFileSync(srcPath, destPath);
                    }
                }
            };
            copyRecursive(publicDir, distPublicDir);
            console.log('📁 Copied public assets to dist/public/');
        }

        // Copy styles directory
        const stylesDir = path.join(__dirname, 'styles');
        if (fs.existsSync(stylesDir)) {
            const copyRecursive = (src, dest) => {
                if (!fs.existsSync(dest)) fs.mkdirSync(dest, { recursive: true });
                for (const entry of fs.readdirSync(src)) {
                    const srcPath = path.join(src, entry);
                    const destPath = path.join(dest, entry);
                    if (fs.statSync(srcPath).isDirectory()) {
                        copyRecursive(srcPath, destPath);
                    } else {
                        fs.copyFileSync(srcPath, destPath);
                    }
                }
            };
            copyRecursive(stylesDir, distDir);
            console.log('📁 Copied styles to dist/');
        }

        // Copy manifest.json for component reference / web app manifest
        const manifestSrc = path.join(__dirname, 'manifest.json');
        if (fs.existsSync(manifestSrc)) {
            fs.copyFileSync(manifestSrc, path.join(distDir, 'manifest.json'));
            console.log('📁 Copied manifest.json to dist/');
        }

        // Copy robots.txt to dist root for SEO
        const robotsSrc = path.join(__dirname, 'public', 'robots.txt');
        if (fs.existsSync(robotsSrc)) {
            fs.copyFileSync(robotsSrc, path.join(distDir, 'robots.txt'));
            console.log('📁 Copied robots.txt to dist/');
        }

        // Copy sitemap.xml to dist root for SEO
        const sitemapSrc = path.join(__dirname, 'sitemap.xml');
        if (fs.existsSync(sitemapSrc)) {
            fs.copyFileSync(sitemapSrc, path.join(distDir, 'sitemap.xml'));
            console.log('📁 Copied sitemap.xml to dist/');
        }

        // Copy .nojekyll to dist root for GitHub Pages
        const nojekyllSrc = path.join(__dirname, '.nojekyll');
        if (fs.existsSync(nojekyllSrc)) {
            fs.copyFileSync(nojekyllSrc, path.join(distDir, '.nojekyll'));
            console.log('📁 Copied .nojekyll to dist/');
        }

        console.log('\n🎉 Build complete!');
    } catch (error) {
        console.error('❌ Build failed:', error.message);
        process.exit(1);
    }
}

buildStaticSite();
