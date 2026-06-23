#!/usr/bin/env python3
"""
Query별 GO 그래프 시각화
E3PQQ8 query의 target GO terms를 색칠해서 표시
"""

import json

def parse_gog_file(filepath):
    """GOG 파일 파싱 - GO term의 hierarchical 정보 저장"""
    go_info = {}  # {go_id: {'name': ..., 'parents': [...], 'category': ...}}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # 형식: GO:ID,description,parent_GO,category
            parts = line.split(',')
            if len(parts) < 3:
                continue
            
            go_id = parts[0]
            description = parts[1]
            parent_info = parts[2] if len(parts) > 2 else ""
            category = parts[3] if len(parts) > 3 else "unknown"
            
            # 부모 정보 추출 (부모는 세미콜론으로 구분될 수 있음)
            parents = []
            if parent_info:
                parent_ids = [p.strip() for p in parent_info.split(';')]
                for p_id in parent_ids:
                    if p_id.startswith('GO:'):
                        parents.append(p_id)
            
            go_info[go_id] = {
                'name': description,
                'parents': parents,
                'category': category
            }
    
    return go_info

def build_subgraph(query_id, target_gos, go_info):
    """Target GO terms의 complete Gene Ontology hierarchy 구성"""
    nodes = {}
    links = []
    
    # Target GO 색상 - orange
    TARGET_GO_COLOR = '#FF8C42'
    
    # 실제 GO root IDs (네임스페이스별)
    # 부모가 없는 GO들이 root입니다
    go_roots = {
        'biological_process': 'GO:0008150',
        'molecular_function': 'GO:0003674',
        'cellular_component': 'GO:0005575'
    }
    
    for ns, go_id in go_roots.items():
        if go_id in go_info:
            go_data = go_info[go_id]
            nodes[go_id] = {
                'id': go_id,
                'label': go_data['name'],
                'title': f"{go_id}: {go_data['name']}",
                'type': 'root',
                'color': '#E0E0E0',
                'category': ns
            }
        else:
            nodes[go_id] = {
                'id': go_id,
                'label': ns.replace('_', ' ').title(),
                'title': f"Root: {ns}",
                'type': 'root',
                'color': '#E0E0E0',
                'category': ns
            }
    
    print(f"  Found GO info for {len(go_info)} terms")
    print(f"  Processing target GOs: {target_gos}")
    print(f"  Creating hierarchical paths...")
    
    target_set = set(target_gos)
    
    # Target GO들에서 root까지 모든 경로 구성
    # visited는 global이 아니라 경로 추적용으로만 사용
    for target_go in target_gos:
        if target_go not in go_info:
            print(f"  ⚠️  GO {target_go} not found in go_info")
            continue
        
        # 각 target마다 독립적으로 경로 추적 (visited를 target별로 초기화)
        visited_in_path = set()
        path = []
        current = target_go
        
        while current:
            if current in visited_in_path:
                break
            
            # Root GO에 도달했으면 경로에 추가하지 않고 링크만 생성
            if current in go_roots.values():
                # 마지막 비-root 노드에서 root GO로 직접 연결
                if path:
                    links.append({
                        'source': path[-1],
                        'target': current,
                        'type': 'subclass_of'
                    })
                break
            
            path.append(current)
            visited_in_path.add(current)
            
            if current not in go_info:
                break
            
            go_data = go_info[current]
            category = go_data.get('category', 'unknown')
            
            # 노드 추가 (중복 가능, 수정 가능)
            is_target = current in target_set
            if is_target:
                color = TARGET_GO_COLOR
                node_type = 'target'
            else:
                color = '#E0E0E0'
                node_type = 'intermediate'
            
            nodes[current] = {
                'id': current,
                'label': go_data['name'],
                'title': f"{current}: {go_data['name']}",
                'type': node_type,
                'color': color,
                'category': category
            }
            
            # 부모로 이동
            parents = go_data.get('parents', [])
            if parents:
                parent_id = parents[0]
                current = parent_id
            else:
                # 부모가 없으면 category(namespace)의 root GO로 연결
                if category in go_roots:
                    root_go_id = go_roots[category]
                    links.append({
                        'source': path[-1],
                        'target': root_go_id,
                        'type': 'subclass_of'
                    })
                current = None
        
        # 경로에 따른 링크 추가
        for i in range(len(path) - 1):
            if path[i] in nodes and path[i+1] in nodes:
                links.append({
                    'source': path[i],
                    'target': path[i+1],
                    'type': 'subclass_of'
                })
    
    return nodes, links

def create_html(query_id, target_gos, nodes, links, output_file):
    """Cytoscape.js를 사용한 인터랙티브 HTML 생성"""
    
    # Cytoscape.js 형식으로 데이터 변환
    cy_elements = []
    
    for node in nodes.values():
        cy_elements.append({
            'data': {
                'id': node['id'],
                'label': node['label'],
                'title': node['title'],
                'color': node['color'],
                'type': node['type']
            }
        })
    
    for link in links:
        cy_elements.append({
            'data': {
                'source': link['source'],
                'target': link['target'],
                'id': f"{link['source']}-{link['target']}",
                'type': link['type']
            }
        })
    
    # 통계 계산
    target_count = len([n for n in nodes.values() if n['type'] == 'target'])
    parent_count = len([n for n in nodes.values() if n['type'] == 'intermediate'])
    root_count = len([n for n in nodes.values() if n['type'] == 'root'])
    total_count = len(nodes)
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Query GO Graph - {query_id}</title>
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
            font-size: 14px;
        }}
        #cy {{
            width: 100%;
            height: 650px;
            background: #f7f7f7;
        }}
        .legend {{
            padding: 20px 30px;
            background: #f9f9f9;
            display: flex;
            gap: 30px;
            flex-wrap: wrap;
            border-bottom: 1px solid #ddd;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .legend-color {{
            width: 25px;
            height: 25px;
            border-radius: 50%;
            border: 2px solid #ddd;
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
            font-family: monospace;
            font-size: 13px;
            line-height: 1.6;
        }}
        .stats {{
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 15px;
            margin-top: 15px;
        }}
        .stat {{
            background: white;
            padding: 12px;
            border-radius: 5px;
            text-align: center;
            border: 1px solid #e0e0e0;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: #667eea;
        }}
        .stat-label {{
            font-size: 12px;
            color: #636e72;
            margin-top: 5px;
        }}
        .controls {{
            padding: 20px 30px;
            background: white;
            display: flex;
            gap: 10px;
            border-bottom: 1px solid #ddd;
            align-items: center;
        }}
        .control-label {{
            font-weight: 600;
            color: #333;
            margin-right: 10px;
        }}
        .toggle-btn {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 5px;
            cursor: pointer;
            font-size: 14px;
            font-weight: 600;
            transition: all 0.3s ease;
            box-shadow: 0 2px 10px rgba(102, 126, 234, 0.3);
        }}
        .toggle-btn:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.5);
        }}
        .toggle-btn:active {{
            transform: translateY(0);
        }}
        .toggle-btn.hidden {{
            background: linear-gradient(135deg, #999 0%, #666 100%);
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Gene Ontology Hierarchy</h1>
            <p>Target GO Terms from Query: <strong>{query_id}</strong></p>
        </div>
        
        <div id="cy"></div>
        
        <div class="controls">
            <span class="control-label">🎮 Visualization Controls:</span>
            <button class="toggle-btn" id="toggle-intermediate-btn" onclick="toggleIntermediateNodes()">
                ✕ Hide Intermediate Nodes
            </button>
        </div>
        
        <div class="legend">
            <div class="legend-item">
                <div class="legend-color" style="background: #FF8C42;"></div>
                <span><strong>Target GO Terms</strong> (E3PQQ8 alignment에서 추출)</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #E0E0E0;"></div>
                <span><strong>Hierarchy Structure</strong> (Root + Intermediate ancestors)</span>
            </div>
            <div class="legend-item">
                <span style="font-size: 12px;">→ subclass_of relationship</span>
            </div>
        </div>
        
        <div class="info">
            <h3>📊 Information</h3>
            <p><strong style="color: #FF8C42;">🟠 Orange nodes</strong> = Target GO terms from E3PQQ8 alignment<br>
               <strong style="color: #999;">⚫ Gray nodes</strong> = Hierarchy structure (root + intermediate ancestors)</p>
            
            <div class="stats">
                <div class="stat">
                    <div class="stat-value">{target_count}</div>
                    <div class="stat-label">Target GO Terms</div>
                </div>
                <div class="stat">
                    <div class="stat-value">{len([n for n in nodes.values() if n['type'] == 'intermediate'])}</div>
                    <div class="stat-label">Intermediate Terms</div>
                </div>
                <div class="stat">
                    <div class="stat-value">{total_count}</div>
                    <div class="stat-label">Total Nodes</div>
                </div>
            </div>
            
            <div class="node-info">
                <strong>Selected Node:</strong> <span id="selected-node">None</span><br>
                <strong>Type:</strong> <span id="node-type">-</span><br>
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
                        'font-size': function(ele) {{
                            return ele.data('type') === 'root' ? '11px' : '10px';
                        }},
                        'width': function(ele) {{
                            if (ele.data('type') === 'root') return '100px';
                            if (ele.data('type') === 'target') return '70px';
                            return '60px';
                        }},
                        'height': function(ele) {{
                            if (ele.data('type') === 'root') return '100px';
                            if (ele.data('type') === 'target') return '70px';
                            return '60px';
                        }},
                        'border-width': function(ele) {{
                            return ele.data('type') === 'target' ? '3px' : '1px';
                        }},
                        'border-color': function(ele) {{
                            return ele.data('type') === 'target' ? '#333' : '#999';
                        }},
                        'text-wrap': 'wrap',
                        'text-max-width': function(ele) {{
                            return ele.data('type') === 'root' ? '90px' : '55px';
                        }},
                        'line-height': '1.2'
                    }}
                }},
                {{
                    selector: 'node:hover',
                    style: {{
                        'background-color': '#ffaa00',
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
                        'line-color': '#999',
                        'target-arrow-color': '#999',
                        'curve-style': 'straight',
                        'width': '2px',
                        'opacity': 0.8
                    }}
                }}
            ],
            layout: {{
                name: 'breadthfirst',
                directed: true,
                roots: ['GO:0008150', 'GO:0003674', 'GO:0005575'],
                spacingFactor: 1.5,
                animate: true,
                animationDuration: 500
            }}
        }});

        // 노드 선택 시 정보 표시
        cy.on('tap', 'node', function(evt) {{
            const node = evt.target;
            document.getElementById('selected-node').textContent = node.data('id');
            document.getElementById('node-type').textContent = node.data('type');
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
                document.getElementById('node-type').textContent = '-';
                document.getElementById('node-description').textContent = '-';
            }}
        }});
        
        // 자동 레이아웃
        setTimeout(() => cy.layout({{ 
            name: 'breadthfirst', 
            directed: true,
            roots: ['GO:0008150', 'GO:0003674', 'GO:0005575']
        }}).run(), 100);
        
        // 중간 노드 보이기/숨기기 토글
        let intermediateNodesVisible = true;
        
        function toggleIntermediateNodes() {{
            const btn = document.getElementById('toggle-intermediate-btn');
            const intermediateNodes = cy.$('node[type="intermediate"]');
            
            if (intermediateNodesVisible) {{
                // 숨기기: intermediate 노드를 통과하는 edge 연결 유지
                
                // 모든 intermediate 노드 처리
                intermediateNodes.forEach(intNode => {{
                    // 이 노드로 들어오는 edge들 (incoming)
                    const incomingEdges = intNode.incomers('edges');
                    // 이 노드에서 나가는 edge들 (outgoing)
                    const outgoingEdges = intNode.outgoers('edges');
                    
                    // incoming source를 outgoing target과 직접 연결
                    incomingEdges.forEach(inEdge => {{
                        const inSource = inEdge.source();
                        outgoingEdges.forEach(outEdge => {{
                            const outTarget = outEdge.target();
                            
                            // 새로운 bypass edge 생성
                            cy.add({{
                                data: {{
                                    id: `bypass-${{inSource.id()}}-${{outTarget.id()}}`,
                                    source: inSource.id(),
                                    target: outTarget.id(),
                                    type: 'bypass'
                                }}
                            }});
                        }});
                    }});
                }});
                
                // Intermediate 노드와 그 직접 연결 edge 숨기기
                intermediateNodes.style('display', 'none');
                cy.elements('edge').forEach(edge => {{
                    const src = edge.source();
                    const tgt = edge.target();
                    if (src.data('type') === 'intermediate' || tgt.data('type') === 'intermediate') {{
                        edge.style('display', 'none');
                    }}
                }});
                
                btn.textContent = '✔ Show Intermediate Nodes';
                btn.classList.add('hidden');
                intermediateNodesVisible = false;
            }} else {{
                // 보이기: bypass edge 제거하고 원래대로 복원
                cy.elements('edge[type="bypass"]').remove();
                cy.elements('edge').style('display', 'element');
                intermediateNodes.style('display', 'element');
                
                btn.textContent = '✕ Hide Intermediate Nodes';
                btn.classList.remove('hidden');
                intermediateNodesVisible = true;
            }}
            
            // 레이아웃 재계산
            setTimeout(() => {{
                const layout = cy.layout({{ 
                    name: 'breadthfirst', 
                    directed: true,
                    roots: ['GO:0008150', 'GO:0003674', 'GO:0005575'],
                    animate: true,
                    animationDuration: 300
                }});
                layout.run();
            }}, 50);
        }}
    </script>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"✅ Gene Ontology Graph HTML 생성 완료: {output_file}")
    print(f"   - Query: {query_id}")
    print(f"   - Target GO terms: {target_count}")
    print(f"   - Intermediate terms: {parent_count}")
    print(f"   - Root namespaces: {root_count}")
    print(f"   - 총 노드: {total_count}")
    print(f"   - 링크: {len(links)}")
    
    # 디버그: nodes 타입별 분류
    print(f"\n  📋 노드 타입별 분류:")
    for node_type in ['root', 'target', 'intermediate']:
        count = len([n for n in nodes.values() if n['type'] == node_type])
        if count > 0:
            print(f"     {node_type}: {count}")

if __name__ == '__main__':
    # E3PQQ8 데이터
    query_id = "E3PQQ8"
    target_gos = ["GO:0090729", "GO:0044556", "GO:0005615", "GO:0005576"]
    
    gog_file = '/Users/yeojingi/Documents/Github/m2g/examples/mmseqs2go/database/uniprot_sprot_func_gog'
    output_file = '/Users/yeojingi/Documents/Github/m2g/examples/mmseqs2go/E3PQQ8_go_graph.html'
    
    print("🔍 GOG 파일 파싱 중...")
    go_info = parse_gog_file(gog_file)
    
    print(f"🎯 Query {query_id}의 subgraph 구성 중...")
    nodes, links = build_subgraph(query_id, target_gos, go_info)
    
    print("🎨 HTML 생성 중...")
    create_html(query_id, target_gos, nodes, links, output_file)
