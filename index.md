---
layout: default
title: Home
---

# Greg Hamilton

**Senior Data Scientist | Ph.D. in Physics | Data Analytics Expert**

Welcome! I'm a data science leader with a unique background bridging theoretical physics, management consulting, and industrial data science. I specialize in applying advanced analytics and machine learning to complex problems in bioinformatics, scientific computing, and business intelligence.

---

## What I Do

🔬 **Bioinformatics & Computational Biology**  
Single-cell RNA sequencing, genomics analysis, and biological data interpretation

🤖 **Machine Learning & AI**  
Model development, deployment, and scaling for real-world applications

📊 **Topological Data Analysis**  
Applying persistent homology and topological methods to complex data structures

🏗️ **Data Science Leadership**  
Building teams, defining strategy, and delivering measurable business impact

---

## Featured Projects

{% assign featured = site.projects | where: "featured", true | sort: "date" | reverse | limit: 3 %}

{% if featured.size > 0 %}
<div class="featured-projects">
  {% for project in featured %}
  <div class="project-card">
    <h3><a href="{{ project.url | relative_url }}">{{ project.title }}</a></h3>
    <p>{{ project.description }}</p>
    <p class="project-tech">
      {% for tech in project.technologies limit:4 %}
        <span class="tech-badge">{{ tech }}</span>
      {% endfor %}
    </p>
    <a href="{{ project.url | relative_url }}" class="button">View Project →</a>
  </div>
  {% endfor %}
</div>

<p style="text-align: center; margin-top: 2em;">
  <a href="{{ '/pages/projects' | relative_url }}" class="button-outline">View All Projects</a>
</p>
{% else %}
<p><em>Featured projects coming soon...</em></p>
{% endif %}

---

## Recent Posts

{% assign recent_posts = site.posts | sort: "date" | reverse | limit: 3 %}

{% if recent_posts.size > 0 %}
<div class="recent-posts">
  {% for post in recent_posts %}
  <div class="post-preview">
    <h3><a href="{{ post.url | relative_url }}">{{ post.title }}</a></h3>
    <p class="post-date">{{ post.date | date: "%B %d, %Y" }}</p>
    {% if post.excerpt %}
    <p>{{ post.excerpt | strip_html | truncatewords: 30 }}</p>
    {% endif %}
  </div>
  {% endfor %}
</div>

<p style="text-align: center; margin-top: 2em;">
  <a href="{{ '/blog' | relative_url }}" class="button-outline">Read More Posts</a>
</p>
{% else %}
<p><em>Blog posts coming soon. Check back for technical tutorials and project deep-dives!</em></p>
{% endif %}

---

## Quick Links

<div class="quick-links">
  <div class="link-section">
    <h3>📋 About</h3>
    <p>Learn about my background, experience, and technical expertise.</p>
    <a href="{{ '/pages/about' | relative_url }}">View Profile →</a>
  </div>

  <div class="link-section">
    <h3>🔬 Projects</h3>
    <p>Explore my portfolio of data science and research projects.</p>
    <a href="{{ '/pages/projects' | relative_url }}">See Projects →</a>
  </div>

  <div class="link-section">
    <h3>⚙️ Pipelines</h3>
    <p>Technical documentation for analysis workflows and pipelines.</p>
    <a href="{{ '/pages/pipelines' | relative_url }}">View Pipelines →</a>
  </div>

  <div class="link-section">
    <h3>📚 Publications</h3>
    <p>Academic papers, technical reports, and presentations.</p>
    <a href="{{ '/pages/publications' | relative_url }}">Read Papers →</a>
  </div>
</div>

---

## Current Focus

I'm currently working on:
- 🧬 Advanced single-cell RNA sequencing analysis methodologies
- 🤖 Machine learning applications in industrial settings
- 📊 Topological data analysis for complex systems
- 👥 Building and mentoring high-performing data science teams

---

## Let's Connect

Interested in collaboration, speaking opportunities, or just want to discuss data science?

<p style="text-align: center; margin-top: 2em;">
  <a href="{{ '/pages/contact' | relative_url }}" class="button">Get in Touch</a>
</p>

---

<div style="text-align: center; color: #888; margin-top: 3em;">
  <p>
    <a href="https://github.com/hreggramDS">GitHub</a> | 
    <a href="#">LinkedIn</a> | 
    <a href="#">Google Scholar</a>
  </p>
</div>
