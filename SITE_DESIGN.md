# Site Design: Greg Hamilton's Data Science Portfolio

## Overview

This site serves as a professional portfolio and research hub for a data science leader, showcasing projects, publications, pipelines, and team deliverables. The site uses Jekyll (GitHub Pages) to provide a clean, maintainable structure that can scale from individual work to team-based research documentation.

## Design Philosophy

**Primary Goals:**
1. **Professional Portfolio** - Showcase technical expertise and project outcomes
2. **Research Documentation** - Detailed write-ups of methodologies, pipelines, and results
3. **Team Collaboration Hub** - Space for documenting team projects with proper attribution
4. **Accessibility** - Easy navigation for recruiters, collaborators, and technical peers

**Target Audiences:**
- Recruiters and hiring managers seeking data science leadership
- Academic collaborators and researchers
- Industry peers interested in methodologies
- Potential team members or collaborators

## Jekyll Configuration Strategy

### Site Structure (Jekyll Conventions)

```
hreggramDS.github.io/
├── _config.yml                 # Site configuration
├── _data/                      # Structured data (YAML/JSON)
│   ├── navigation.yml          # Site navigation structure
│   ├── projects.yml            # Project metadata
│   └── publications.yml        # Publications list
├── _includes/                  # Reusable components
│   ├── header.html             # Site header
│   ├── footer.html             # Site footer
│   ├── project-card.html       # Project display component
│   └── figure.html             # Figure with caption component
├── _layouts/                   # Page templates
│   ├── default.html            # Base layout
│   ├── page.html               # Standard page
│   ├── post.html               # Blog post layout
│   └── project.html            # Project showcase layout
├── _posts/                     # Blog posts (date-based)
│   └── YYYY-MM-DD-title.md     # Jekyll naming convention
├── _projects/                  # Project collection (custom)
│   └── project-name.md         # Individual projects
├── _pipelines/                 # Pipeline documentation (custom)
│   └── pipeline-name.md        # Pipeline descriptions
├── assets/                     # Static assets
│   ├── css/                    # Stylesheets
│   ├── js/                     # JavaScript
│   ├── images/                 # Images
│   │   ├── projects/           # Project-specific images
│   │   ├── figures/            # Research figures
│   │   └── profile/            # Profile photos
│   └── files/                  # Downloadable files (PDFs, etc.)
├── pages/                      # Static pages
│   ├── about.md                # About/CV page
│   ├── projects.md             # Projects index
│   ├── publications.md         # Publications list
│   ├── pipelines.md            # Pipelines index
│   └── contact.md              # Contact information
├── index.md                    # Landing page
└── README.md                   # Repository documentation
```

## Content Organization

### Collections Strategy

Define custom collections in `_config.yml` for different content types:

```yaml
collections:
  projects:
    output: true
    permalink: /projects/:name/
  pipelines:
    output: true
    permalink: /pipelines/:name/
  team:
    output: true
    permalink: /team/:name/
```

### Page Types & Purpose

#### 1. **Landing Page** (`index.md`)
- Brief introduction and value proposition
- Featured projects (3-4 cards)
- Quick links to major sections
- Recent posts/updates
- Professional headshot and tagline

#### 2. **About/CV Page** (`pages/about.md`)
- Professional biography (expanded from current about.md)
- Education and credentials
- Work experience timeline
- Technical skills matrix
- Downloadable CV/resume
- Contact information

#### 3. **Projects Page** (`pages/projects.md`)
- Grid/card layout of all projects
- Filterable by: domain, technology, status, team size
- Each project card includes:
  - Title and brief description
  - Key technologies
  - Visual thumbnail
  - Link to detailed project page

#### 4. **Individual Project Pages** (`_projects/*.md`)
Each project should follow this structure:
- **Header**: Title, date, team members, status
- **Overview**: Problem statement and objectives
- **Methodology**: Technical approach and architecture
- **Implementation**: Key code snippets, tools, technologies
- **Results**: Metrics, visualizations, outcomes
- **Figures**: Charts, diagrams, screenshots with captions
- **Deliverables**: Links to repos, papers, dashboards
- **Learnings**: Insights and takeaways

#### 5. **Pipelines Page** (`pages/pipelines.md`)
- Documentation hub for data pipelines and workflows
- Organized by domain (e.g., scRNA analysis, NLP, computer vision)
- Each pipeline includes:
  - Architecture diagrams
  - Step-by-step process flow
  - Code examples and snippets
  - Performance benchmarks
  - Usage instructions
  - Links to implementation repos

#### 6. **Blog Posts** (`_posts/`)
- Technical write-ups and tutorials
- Project deep-dives
- Research updates
- Conference notes
- Methodology discussions
- Use Jekyll's date-based naming: `YYYY-MM-DD-title.md`

#### 7. **Publications Page** (`pages/publications.md`)
- Academic papers and articles
- Conference presentations
- Technical reports
- Organized chronologically with:
  - Title and authors
  - Publication venue and date
  - Abstract/summary
  - Links to PDF, DOI, repository

## Navigation Structure

### Primary Navigation
- **Home** - Landing page
- **About** - Professional background and CV
- **Projects** - Portfolio of work
- **Pipelines** - Technical documentation
- **Publications** - Research outputs
- **Blog** - Technical posts
- **Contact** - Get in touch

### Secondary Navigation (Footer)
- Social links (GitHub, LinkedIn, Google Scholar)
- RSS feed for blog
- Site map
- Privacy/disclaimers if needed

## Content Guidelines

### Front Matter Templates

#### Project Template
```yaml
---
layout: project
title: "Project Title"
date: 2025-06-15
team: ["Greg Hamilton", "Collaborator Name"]
status: "Completed" # or "In Progress", "Published"
domain: ["Machine Learning", "NLP"]
technologies: ["Python", "PyTorch", "Transformers"]
github: "https://github.com/username/repo"
thumbnail: "/assets/images/projects/project-name.png"
description: "Brief one-line description"
featured: true # Show on landing page
---
```

#### Pipeline Template
```yaml
---
layout: pipeline
title: "scRNA-seq Analysis Pipeline"
date: 2025-01-20
version: "2.0"
domain: "Bioinformatics"
technologies: ["Python", "Scanpy", "scanpy-workflow"]
repository: "https://github.com/username/scrna-pipeline"
documentation: "https://docs.example.com"
---
```

#### Post Template
```yaml
---
layout: post
title: "Understanding Many-Body Localization"
date: 2025-02-15
categories: ["physics", "research"]
tags: ["quantum", "topology", "data-analysis"]
author: "Greg Hamilton"
excerpt: "A brief excerpt for previews"
---
```

### Writing Style Guidelines

**For Technical Audiences:**
- Lead with the problem and solution
- Include code snippets with syntax highlighting
- Provide architecture diagrams and flowcharts
- Link to GitHub repositories and live demos
- Use mathematical notation (LaTeX) when appropriate

**For General Audiences:**
- Start with business impact and outcomes
- Use visualizations over equations
- Explain technical concepts in plain language
- Highlight metrics and measurable results

**Visual Elements:**
- Use high-quality figures (minimum 1200px width for main images)
- Include alt text for accessibility
- Provide captions with figure numbers for reference
- Create consistent visual style (color palette, fonts)

## Jekyll Features to Leverage

### 1. **Liquid Templating**
- Use loops to generate project grids
- Conditional display based on `featured` flag
- Dynamic sidebar with recent posts
- Auto-generate table of contents

### 2. **Data Files** (`_data/`)
- Store structured data separately (projects, publications)
- Easier to update than markdown front matter
- Enable consistent formatting across site

### 3. **Includes for Reusable Components**
```liquid
{% include project-card.html 
   title="Project Name"
   image="/assets/images/project.png"
   description="Description text"
   link="/projects/project-name/" %}
```

### 4. **Syntax Highlighting**
- Configure Rouge or Pygments for code blocks
- Support multiple languages (Python, R, SQL, bash)
- Add line numbers for longer snippets

### 5. **Pagination**
- Paginate blog posts for better performance
- Configure posts per page (10-15 recommended)

### 6. **Tags and Categories**
- Tag posts by topic/technology
- Create tag archive pages
- Enable filtering and discovery

## Design and Branding

### Visual Theme
- **Current**: Jekyll Theme Slate (dark, professional)
- **Consider**: Custom theme building on Slate
- **Key Elements**:
  - Clean, readable typography
  - Data visualization-friendly color palette
  - Professional but approachable tone
  - Mobile-responsive design

### Color Palette (Data Science Theme)
- Primary: Deep blue (#2C3E50) - trust, intelligence
- Accent: Bright cyan (#3498DB) - data, technology
- Highlight: Orange (#E74C3C) - action, energy
- Neutral: Grays for text and backgrounds
- Success: Green for metrics and positive results

### Typography
- Headers: Clean sans-serif (Source Sans Pro, Inter, Roboto)
- Body: Readable serif or sans-serif (Georgia, Lora, Open Sans)
- Code: Monospace (Fira Code, Source Code Pro)

## Implementation Phases

### Phase 1: Foundation (Week 1)
- [x] Basic _config.yml setup
- [ ] Create proper directory structure
- [ ] Set up collections (projects, pipelines)
- [ ] Migrate existing content to new structure
- [ ] Create basic layouts (default, page, post, project)

### Phase 2: Core Content (Week 2-3)
- [ ] Develop comprehensive About/CV page
- [ ] Create project showcase pages
- [ ] Document key pipelines (starting with scRNA)
- [ ] Set up navigation and site structure
- [ ] Import or create 3-5 featured projects

### Phase 3: Enhancement (Week 4)
- [ ] Add data files for projects and publications
- [ ] Create reusable includes (project cards, figure components)
- [ ] Implement filtering and tagging system
- [ ] Add social sharing capabilities (Greg's note: omit this)
- [ ] Set up Google Analytics or privacy-friendly alternative

### Phase 4: Content Population (Ongoing)
- [ ] Migrate all historical projects
- [ ] Write blog posts on key topics
- [ ] Document all major pipelines
- [ ] Add publication list
- [ ] Regular updates and maintenance

## Jekyll-Specific Best Practices

### Referencing Content

**Posts:**
```markdown
See my [previous post]({% post_url 2025-02-15-scrna-analysis %}) for details.
```

**Pages:**
```markdown
Check out [my projects]({{ site.baseurl }}{% link pages/projects.md %})
```

**Assets:**
```markdown
![Architecture Diagram]({{ site.baseurl }}/assets/images/projects/architecture.png)
```

**Collections:**
```liquid
{% for project in site.projects %}
  <div>{{ project.title }}</div>
{% endfor %}
```

### Performance Optimization
- Compress images (use WebP when possible)
- Minimize CSS and JavaScript
- Use lazy loading for images
- Enable browser caching
- Consider CDN for assets

### SEO and Discoverability
- Add meta descriptions to all pages
- Use semantic HTML5 structure
- Include Open Graph tags for social sharing
- Submit sitemap to search engines
- Use descriptive URLs (permalinks)

## Maintenance and Growth

### Content Calendar
- **Weekly**: New blog post or project update
- **Monthly**: Pipeline documentation or tutorial
- **Quarterly**: Update CV and project portfolio
- **Annually**: Site design refresh and audit

### Analytics to Track
- Most visited projects
- Popular blog posts
- Traffic sources
- User journey and navigation patterns
- Bounce rate and engagement

### Future Enhancements
- Interactive visualizations (D3.js, Plotly)
- Jupyter notebook integration
- Search functionality
- Comments system (if desired)
- Newsletter signup
- Team member profiles section
- Client testimonials (if applicable)

## Technical Stack Summary

**Core:**
- Jekyll 4.x (GitHub Pages compatible)
- Liquid templating language
- Markdown for content
- YAML for configuration

**Styling:**
- CSS/SCSS (Sass support built into Jekyll)
- Responsive design (mobile-first)
- Custom overrides to Jekyll Theme Slate

**JavaScript (Minimal):**
- Navigation toggle for mobile
- Project filtering
- Image lightbox (optional)
- Syntax highlighting

**Hosting:**
- GitHub Pages (free, automatic deployment)
- Custom domain (optional): greghamiton.com or similar

## Migration Plan from Current Structure

### Immediate Actions:
1. Move `contact/about.md` → `pages/about.md`
2. Rename `posts/` → `_posts/` and update filenames to Jekyll convention
3. Create `_projects/` and migrate `scRNA_paper.md` as first project
4. Update `_config.yml` with collections and site metadata
5. Create basic layouts in `_layouts/`
6. Update README.md to focus on site usage (move test content to actual pages)

### Content Transformation:
- **Current** `posts/scRNA_paper.md` → **New** `_projects/scrna-analysis.md` (full project page)
- **Current** `contact/about.md` → **New** `pages/about.md` (enhanced CV page)
- **Current** README test content → **New** `index.md` (professional landing page)

---

## Next Steps

1. Review this design document
2. Approve overall structure and approach
3. Begin Phase 1 implementation:
   - Create directory structure
   - Set up collections
   - Create basic layouts
4. Prioritize which projects to showcase first
5. Gather assets (images, figures, papers) for content population

This design provides a scalable, professional foundation that can grow from a personal portfolio to a team research hub while maintaining Jekyll's simplicity and GitHub Pages' hosting benefits.
