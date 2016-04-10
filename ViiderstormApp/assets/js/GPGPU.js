//TODO:
// -add change listener to number of scene elements and recompile shaders if it changes
// -add change listener to scene object origins so objects can be moved around without recompiation


function GPGPU(SCR, width, height, rt){
    var self = this;
    
    //fields
    self.rayTracer          = rt;
    self.intersectionShader = null;
    self.colorShader        = null;
    self.objectShaderCode   = [];
    self.gpRenderer         = SCR.scRenderer;
    self.gpCamera           = (function(){
        var i      = new THREE.OrthographicCamera( 1 / -2, 1 / 2, 1 / 2, 1 / -2, 1, 1000 )
        i.position = new THREE.Vector3(0,0,1);
        return i;
    })()
    self.gpScene             = (function(){
        var sc       = new THREE.Scene();
        var planeGeo = new THREE.PlaneGeometry(1,1);
        var plane    = new THREE.Mesh(planeGeo, new THREE.MeshBasicMaterial());
        plane.position.z = -1;
        //plane.rotation.z = -1.5708;
        sc.add(plane);
        return sc;
    })()
    
    //init
    self.init                = function(){
        var tmp = [];
        for(var i = 0; i < self.rayTracer.scene.geometries.length; i++){
            var s = self.rayTracer.scene.geometries[i].getShaderCode(i);
            tmp.push(s);
        }
        self.objectShaderCode   = tmp;
        self.intersectionShader = self.getIntersectionShader(self.rayTracer.scene);
        self.colorShader        = self.getColorShader(self.rayTracer.scene);
    }
    
    //class functions
    self.getIntersections = function(rayData){
        self.updateUniforms(self.intersectionShader, rayData );
        var intersections = self.compute(self.intersectionShader); // vec4(x,y,z, objectIndex)
        return intersections;
    }
    
    self.getColors        = function(intersectionData){
        self.updateUniforms(self.colorShader, intersectionData);
        var colors        = self.compute(self.colorShader);
        return colors;
    }
    
    self.compute          = function(shaderMat){
        self.gpScene.children[0].material = shaderMat;
        var output = new THREE.WebGLRenderTarget(width, height, {
            minFilter: THREE.NearestFilter,
			magFilter: THREE.NearestFilter,
			format: THREE.RGBAFormat,
			type: THREE.FloatType
        });
        self.gpRenderer.render(self.gpScene, self.gpCamera, output);
        return output;
    }
    
    self.updateUniforms        = function(shader, uniforms){
        shader.uniforms = $.extend(shader.uniforms, uniforms);
    }
    
    self.getIntersectionShader = function(scene){
        var sceneElems = scene.geometries;
        var shader = "precision mediump float;\nuniform sampler2D origins;\nuniform sampler2D directions;\nvarying vec2 vUv;\n";
        
        for(var i = 0 ; i < self.objectShaderCode.length; i++){
            shader += self.objectShaderCode[i].boilerCode;
        }
        shader += "vec4 findClosestIntersection(vec4[" + self.objectShaderCode.length +"] intersections){\n" +
                  "\tvec4 max = vec4(9999999,9999999,9999999, -1.0);\n" +
                  "\tfor(int i = 0; i < " + self.objectShaderCode.length + "; i++){\n"+
                  "\t\tif(intersections[i].w < 0.0){\n"+
                  "\t\t\tcontinue;\n"+
                  "\t\t}\n"+
                  "\t\tfloat distance1 = distance(intersections[i].xyz, texture2D(origins, vUv.yx).xyz);\n" +
                  "\t\tfloat distance2 = distance(max.xyz, texture2D(origins, vUv.yx).xyz);\n" +
                  "\t\tif(distance1 < distance2){\n"+
                  "\t\t\tmax = intersections[i].xyzw;\n" +
                  "\t\t}\n" +
                  "\t}\n" +
                  "\treturn max;\n" +
                  "}\n"
        for(var i = 0 ; i < self.objectShaderCode.length; i++){
            shader += self.objectShaderCode[i].functionCode;
        }
        shader += "void main(){\n";
        for(var i = 0; i < self.objectShaderCode.length; i++){
            if(sceneElems[i] instanceof Plane){
                shader += self.objectShaderCode[i].faceDefine;
            }
        }
        shader += "\tvec4 intersections[" + self.objectShaderCode.length + "];\n"
        for(var i = 0; i < self.objectShaderCode.length; i++){
            var qual = self.objectShaderCode[i].qualifier;
            if(sceneElems[i] instanceof Plane){
                shader += "\tintersections[" + i + "] = " + qual("getIntersection") + 
                "(texture2D(origins,vUv.yx).xyz,texture2D(directions,vUv.yx).xyz," + qual("faces") + ");\n";
            } else {
                shader += "\tintersections[" + i + "] = " + qual("getIntersection") + "(texture2D(origins,vUv.yx).xyz,texture2D(directions,vUv.yx).xyz);\n";
            }
        }
        
        shader += "\tvec4 closestIntersection = findClosestIntersection(intersections);\n" +
                  "\tif(closestIntersection.w >= 0.0){\n" +
                  "\t\tgl_FragColor = closestIntersection;\n"+
                  "\t} else {\n" +
                  "\t\tgl_FragColor = vec4(0,0,0,9999999);\n"+ //use a exceptionally large number instead of negativies. Negatives don't return correctly.
                  "\t}\n" +
                  "}\n";
        
        var uniforms = {};
        for(var i = 0 ; i < self.objectShaderCode.length; i++){
            uniforms = $.extend(uniforms, self.objectShaderCode[i].intersect_uniforms);
        }
        
        return new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: document.getElementById("passThruShader").textContent,
            fragmentShader: shader
        });
    }
    
    self.getColorShader   = function(scene){
        var sceneElems = scene.geometries;
        var shader = "precision mediump float;\nuniform sampler2D intersections;\nvarying vec2 vUv;\n"; //TODO: Add normal and light info
        shader += "uniform vec3 bgColor;\n"
        for(var i = 0 ; i < self.objectShaderCode.length; i++){
            shader += "uniform vec3 " + self.objectShaderCode[i].qualifier("color") + ";\n";
            shader += "uniform float " + self.objectShaderCode[i].qualifier("id") + ";\n";
        }
        shader += "void main(){\n";
        shader += "\tvec4 intersection = texture2D(intersections, vUv.xy);\n"
        for(var i = 0 ; i < self.objectShaderCode.length; i++){
            var qual = self.objectShaderCode[i].qualifier;
            var begin = (i == 0) ? "\tif" : "else if";
            begin += "( intersection.w == " + qual("id") + "){\n";
            begin += "\t\tgl_FragColor = vec4(" + qual("color") + ", 1.0);\n";
            begin += "\t} ";
            shader += begin;
        }
        shader += "\telse {\n";
        shader += "\t\tgl_FragColor = vec4(bgColor, 1.0);\n";
        shader += "\t}\n";
        shader += "}";
        var uniforms = {}
        for(var i = 0 ; i < self.objectShaderCode.length; i++){
            uniforms = $.extend(uniforms, self.objectShaderCode[i].color_uniforms);
        }
        
        return new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: document.getElementById("passThruShader").textContent,
            fragmentShader: shader
        });
        
    }
    
    self.init();
}