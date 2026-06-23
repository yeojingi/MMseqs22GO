#!/usr/bin/env python3
"""
GO (Gene Ontology) 그래프 시각화 HTML 생성 스크립트
uniprot_sprot_func_gog 파일에서 샘플 데이터를 추출해 상호작용하는 그래프로 표시
"""

import json
import sys

def parse_go_file(filepath, limit=20):
    """GO 파일을 파싱해서 그래프 데이터 생성"""
    nodes = {}
    links = []
    node_id = 0
    
    with open(filepath, 'r') as f:
        for line_idx, line in enumerate(f):
            if line_idx >= limit:
                break
            
            line = line.strip()
            if not line:
                continue
            
            # 형식: GO:ID,description,parent_GO;more_parent,category
            parts = line.split(',')
            if len(parts) < 2:
                continue
            
            go_id = parts[0]
            description = parts[1]
            category = parts[-1] if len(parts) > 2 else "unknown"
            
            # 현재 노드 추가
            if go_id not in nodes:
                nodes[go_id] = {
                    'id': go_id,
                    'label': f"{go_id}\n{description[:30]}...",
                    'title': f"{go_id}: {description}",
                    'category': category,
                    'color': get_color_by_category(category)
                }
            
            # 부모 노드 관계 파싱
            parent_info = parts[2] if len(parts) > 2 else ""
            if parent_info and not parent_info.startswith("GO:"):
                # 마지막이 category이면 그 앞이 parent_info
                if len(parts) > 3:
                    parent_info = parts[2]
            
            if parent_info and parent_info != category:
                parent_gos = [p.strip() for p in parent_info.split(';')]
                for parent_go in parent_gos:
                    if parent_go.startswith('GO:'):
                        if parent_go not in nodes:
                            nodes[parent_go] = {
                                'id': parent_go,
                                'label': parent_go,
                                'title': parent_go,
                                'category': 'parent',
                                'color': '#cccccc'
                            }
                        links.append({
                            'source': parent_go,
                            'target': go_id,
                            'type': 'parent_of'
                        })
    
    return list(nodes.values()), links

def get_color_by_category(category):
    """카테고리별 색상 결정"""
    colors = {
        'biological_process': '#FF6B6B',
        'molecular_function': '#4ECDC4',
        'cellular_component': '#45B7D1',
        'parent': '#cccccc',
        'unknown': '#95a5a6'
    }
    return colors.get(category, '#95a5a6')

def create_html(nodes, links, output_file):
    """Cytoscape.js를 사용한 인터랙티브 HTML 생성"""
    
    # Cytoscape.js 형식으로 데이터 변환
    cy_elements = []
    
    for node in nodes:
        cy_elements.append({
            'data': {
                'id': node['id'],
                'label': node['label'],
                'title': node['title'],
                'color': node['color']
            }
        })
    
    for link in links:
        cy_elements.append({
            'data': {
                'source': link['source'],
                'target': link['target'],
                'id': f"{link['source']}-{link['target']}"
            }
        })
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Gene Ontology Graph Visualization</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.20.1/cytoscape.min.js"></script>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            margin: 0;
            font-size: 28px;
        }}
        .header p {{
            margin: 10px 0 0 0;
            opacity: 0.9;
        }}
        #cy {{
            width: 100%;
            height: 600px;
            background: #f7f7f7;
        }}
        .legend {{
            padding: 20px 30px;
            background: #f9f9f9;
            display: flex;
            gap: 30px;
            flex-wrap: wrap;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .legend-color {{
            width: 20px;
            height: 20px;
            border-radius: 50%;
            border: 1px solid #ddd;
        }}
        .info {{
            padding: 30px;
            background: #f0f4ff;
            border-left: 4px solid #667eea;
        }}
        .info h3 {{
            margin-top: 0;
            color: #667eea;
        }}
        .node-info {{
            background: white;
            padding: 15px;
            margin-top: 10px;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
            min-height: 20px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Gene Ontology Graph Visualization</h1>
            <p>Sample GO terms from Swiss-Prot database with hierarchical relationships</p>
        </div>
        
        <div id="cy"></div>
        
        <div class="legend">
            <div class="legend-item">
                <div class="legend-color" style="background: #FF6B6B;"></div>
                <span>Biological Process</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #4ECDC4;"></div>
                <span>Molecular Function</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #45B7D1;"></div>
                <span>Cellular Component</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #cccccc;"></div>
                <span>Parent Term</span>
            </div>
        </div>
        
        <div class="info">
            <h3>📊 Information</h3>
            <p>Click on nodes to highlight relationships. Drag nodes to rearrange.</p>
            <div class="node-info">
                <strong>Selected Node:</strong> <span id="selected-node">None</span><br>
                <strong>Description:</strong> <span id="node-description">-</span>
            </div>
        </div>
    </div>

    <script>
        const cy = cytoscape({{
            container: document.getElementById('cy'),
            elements: {json.dumps(cy_elements)},
            style: [
                {{
                    selector: 'node',
                    style: {{
                        'content': 'data(label)',
                        'background-color': 'data(color)',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'font-size': '11px',
                        'width': '60px',
                        'height': '60px',
                        'border-width': '2px',
                        'border-color': '#333',
                        'text-wrap': 'wrap',
                        'text-max-width': '55px'
                    }}
                }},
                {{
                    selector: 'node:hover',
                    style: {{
                        'background-color': '#ff0000',
                        'box-shadow': '0 0 20px rgba(0,0,0,0.5)'
                    }}
                }},
                {{
                    selector: 'node:selected',
                    style: {{
                        'border-width': '3px',
                        'border-color': '#000',
                        'box-shadow': '0 0 30px rgba(102, 126, 234, 0.5)'
                    }}
                }},
                {{
                    selector: 'edge',
                    style: {{
                        'target-arrow-shape': 'triangle',
                        'line-color': '#ccc',
                        'target-arrow-color': '#ccc',
                        'curve-style': 'bezier',
                        'width': '2px'
                    }}
                }}
            ],
            layout: {{
                name: 'cose',
                directed: true,
                padding: 10,
                animate: true,
                animationDuration: 500
            }}
        }});

        // 노드 선택 시 정보 표시
        cy.on('tap', 'node', function(evt) {{
            const node = evt.target;
            document.getElementById('selected-node').textContent = node.data('id');
            document.getElementById('node-description').textContent = node.data('title');
            
            // 다른 노드들을 흐리게 처리
            cy.elements().removeClass('faded');
            cy.elements().not(node.outgoers().union(node.incomers()).union(node)).addClass('faded');
        }});

        // 배경 클릭 시 초기화
        cy.on('tap', function(evt) {{
            if (evt.target === cy) {{
                cy.elements().removeClass('faded');
                document.getElementById('selected-node').textContent = 'None';
                document.getElementById('node-description').textContent = '-';
            }}
        }});
    </script>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"✅ HTML 파일 생성 완료: {output_file}")
    print(f"   - 노드 수: {len(nodes)}")
    print(f"   - 링크 수: {len(links)}")

if __name__ == '__main__':
    go_file = '/Users/yeojingi/Documents/Github/m2g/examples/mmseqs2go/database/uniprot_sprot_func_gog'
    output_file = '/Users/yeojingi/Documents/Github/m2g/examples/mmseqs2go/go_graph_sample.html'
    
    print("🔍 GO 파일 파싱 중...")
    nodes, links = parse_go_file(go_file, limit=20)
    
    print("🎨 HTML 생성 중...")
    create_html(nodes, links, output_file)
