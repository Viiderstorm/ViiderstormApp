
/*
* int iwidth - width of imageplane
* int iheight - height of imageplane
* THREE.Vector3 location - origin point of the ImageFrame
*/
function ImageFrame(iwidth, iheight, location) {
    
    /* Instance Variables */
    var self = this;
    self.imageWidth = iwidth;
    self.imageHeight = iheight;

    self.planeWidth = 2;
    self.planeHeight = 2;
    self.location = location; //three vector3

    self.pixels = (function () {
        var tmp = [];
        for (var i = 0; i < self.imageHeight; i++) {
            tmp.push([])
            for (var k = 0; k < self.imageWidth; k++) {
                tmp[i].push(new THREE.Color(0x000000))
            }
        }
        return tmp
    })()

    self.pixelWidth = self.planeWidth / self.imageWidth;
    self.pixelHeight = self.planeHeight / self.imageHeight;

    /* generateImageData
    * Description: Converts ray data into pixel data that the HTML canvas element can understand
    */
    self.generateImageData = function () {
        var tmp = [];
        for (var i = 0; i < self.imageHeight; i++) {
            for (var k = 0; k < self.imageWidth; k++) {

                var style = self.pixels[i][k].getStyle().replace("rgb(", "").replace(")", "")
                var values = style.split(",");

                tmp.push(parseInt(values[0]));
                tmp.push(parseInt(values[1]));
                tmp.push(parseInt(values[2]));
                tmp.push(255);
            }
        }

        return new Uint8ClampedArray(tmp);
    }

    /* displayImage
    * Description: Takes HTML Canvas Image Data and places it in the canvas element
    *
    * sstring canvasId - HTML ID of the canvas to draw the image onto.
    */
    self.displayImage = function (canvasId) {

        var c = document.getElementById(canvasId);
        var ctx = c.getContext("2d");
        var imgData = ctx.createImageData(self.imageHeight, self.imageWidth);
        imgData.data.set(self.generateImageData());
        ctx.putImageData(imgData, 0, 0);

    }

}


/* Class RayTracer
* THREE.Camera camera - THREE js camera object
* ImagePlane   frame  - ImagePlane used to capture the Ray Traced Image
* THREE.Scene  scene  - Scene to use the Ray Tracer on  
*/
function RayTracer(camera, frame, scene) {
    var self = this;

    self.frame = frame;
    self.geometries = (function () {
        var tmp = []
        for (var i = 0; i < scene.children.length; i++) {
            if (scene.children[i] instanceof THREE.Mesh) {
                tmp.push(scene.children[i]);
            }
        }
        return tmp;
    })()

    self.lights = (function () {
        var tmp = []
        for (var i = 0; i < scene.children.length; i++) {
            if (scene.children[i] instanceof THREE.PointLight) {
                tmp.push(scene.children[i]);
            }
        }
        return tmp;
    })()

    self.rays = generateRays(camera, frame);

    self.traceRays = function () {
        for (var i = 0; i < self.rays.length; i++) {
            for (var k = 0; k < self.rays[i].length; k++) {

                var intersectedObject = findIntersection(self.rays[i][k], self.geometries);
                if (intersectedObject != null) {
                    self.frame.pixels[i][k] = intersectedObject.material.color;
                }

            }
        }

        return self.frame;

    }

    self.render = function (canvasId) {
        var imgPlane = self.traceRays();
        imgPlane.displayImage(canvasId);
    }

}

/* Class Ray
* THREE.Vector3 origin    - origin of the ray
* THREE.Vector3 direction - direction of the ray
*/
function Ray(origin, direction) {
    var self = this;

    self.origin = origin;
    self.direction = direction;
}

/* Function generateRays
* Description: Takes a Camera and an ImagePlane and calculates all the ray origins and direction for each pixel of the ImagePlane
*
* THREE.Camera camera     - Camera to use for the ray tracer
* ImagePlane   frame      - ImagePlane that will capture the ray traced image
*/
function generateRays(camera, frame) {
    
    //assume camera is centered behind frame, looking through the center of frame
    var cameraLocation = camera.position.clone();
    var originalFramePixel = frame.location.clone();
    originalFramePixel.x = originalFramePixel.x - (frame.planeWidth / 2) + (0.5 * frame.pixelWidth);
    originalFramePixel.y = originalFramePixel.y + (frame.planeHeight / 2) - (0.5 * frame.pixelHeight);
    var currentFramePixel = originalFramePixel.clone();

    var rays = [];
    for (var i = 0; i < frame.imageHeight; i++) {
        rays.push([]);
        for (var k = 0; k < frame.imageWidth; k++) {

            var rayDirection = currentFramePixel.clone().sub(cameraLocation).normalize();

            rays[i].push(new Ray(currentFramePixel.clone(), rayDirection))

            currentFramePixel.x += frame.pixelWidth;

        }

        currentFramePixel = originalFramePixel.clone();
        currentFramePixel.y -= (frame.pixelHeight * (i + 1));
    }

    return rays;

}

/* Function findIntersections
* Description: Finds the intersection points given a ray and a list of THREE.Geometry
* Ray ray - ray to find intersections for
* Geometry[] geometries - geometries to check against
*/
function findIntersection(ray, geometries) {

    var intersections = [];

    for (var i = 0; i < geometries.length; i++) {
        if (geometries[i].geometry instanceof THREE.SphereGeometry) {
            var intersectPoint = getSphereIntersection(ray, geometries[i]);
            if (intersectPoint != null) {
                intersections.push({
                    geometry: geometries[i],
                    intersection: intersectPoint
                });
            }
        }
        if (geometries[i].geometry instanceof THREE.PlaneGeometry) {
            var intersectPoint = getPlaneIntersectionNotSchool(ray, geometries[i]);
            if (intersectPoint != null && validatePointIsInPolygonPlane(intersectPoint, geometries[i])) {
                intersections.push({
                    geometry: geometries[i],
                    intersection: intersectPoint
                });
            }
        }
    }

    var closestIntersectionDistance = 99999;
    var closestGeometry = null;
    for (var i = 0; i < intersections.length; i++) {

        var distance = intersections[i].intersection.clone().sub(ray.origin).length();
        if (distance < closestIntersectionDistance) {
            closestIntersectionDistance = distance;
            closestGeometry = intersections[i].geometry;
        }

    }

    return closestGeometry;

}

/* Function getSphereIntersection
* Decription: Calculates intersection with sphere, returns null if no intersection is found.
* Ray ray - ray to find an intersection with
* THREE.Sphere - sphere to check for intersection against
*/
function getSphereIntersection(ray, sphere) {

    var dx = ray.direction.x;
    var dy = ray.direction.y;
    var dz = ray.direction.z;

    var xc = sphere.position.x;
    var yc = sphere.position.y;
    var zc = sphere.position.z;

    var xo = ray.origin.x;
    var yo = ray.origin.y;
    var zo = ray.origin.z;

    var r = sphere.geometry.parameters.radius;

    var A = Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2);
    var B = 2 * ((dx * (xo - xc)) + (dy * (yo - yc)) + (dz * (zo - zc)));
    var C = Math.pow(xo - xc, 2) + Math.pow(yo - yc, 2) + Math.pow(zo - zc, 2) - Math.pow(r, 2);

    var W1 = (-B - Math.sqrt(Math.pow(B, 2) - (4 * A * C))) / (2 * A);
    var W2 = (-B + Math.sqrt(Math.pow(B, 2) - (4 * A * C))) / (2 * A);

    var W = W1 > 0 ? W1 : W2;

    if (W >= 0) {
        return ray.origin.clone().add(ray.direction.clone().multiplyScalar(W));
    }

    return null;

}


/*
* Equation for plane intersection from class. Very confusing, didn't end up using it in the end.
*/
function getPlaneIntersection(ray, plane) {

    var pointOnPlace = ray.origin.clone().sub(plane.geometry.vertices[0]);
    var v = plane.geometry.vertices.slice(0,4);
    var normal = getNormal([v[0], v[2], v[1]]).normalize();
    var f = normal.dot(plane.geometry.vertices[0].clone().sub(new THREE.Vector3(0, 0, 0))) / normal.length();
    var denom = normal.dot(ray.direction)
    var W = -((normal.dot(pointOnPlace) + f) / denom)

    if (W > 0) {
        return ray.origin.clone().add(ray.direction.clone().multiplyScalar(W));
    }

    return null;

}


/*
* Equation for plane intersection I found on the web. Made more sense to me and worked very well.
* Ray ray - intersecting ray
* Three.Mesh plane - plane to find intersection with 
*/
function getPlaneIntersectionNotSchool(ray, plane){
    
    var pointOnPlane = plane.geometry.vertices[0];
    var v = plane.geometry.vertices.slice(0,4);
    var normal = getNormal([v[0], v[2], v[1]]).normalize();
    var rayToPoint = ray.origin.clone().sub(pointOnPlane);
    var denom = normal.clone().dot(ray.direction);
    var num = normal.clone().dot(rayToPoint);
    
    var t = -(num / denom);
    if(t > 0){
        return ray.origin.clone().add(ray.direction.clone().multiplyScalar(t));
    }
    return null;
    
}

/*
* Description: Checks that the intersection point lies within one of the faces of the plane
* THREE.Vector3 point - Intersection point
* THREE.Mesh polygon - plane to intersect with
*/
function validatePointIsInPolygonPlane(point, polygon) {

    var successfulFaceTests = 0;
    for (var i = 0; i < polygon.geometry.faces.length; i++) {
        if (checkPointInTriangle(point, polygon.geometry.vertices, polygon.geometry.faces[i])) {
            successfulFaceTests += 1;
        }
    }

    if (successfulFaceTests > 0) {
        return true;
    }

    return false;


}


/*
* Description: Checks that the intersection point lies within a face
* THREE.Vector3 point - Intersection point
* THREE.Vector3[] vertices- vertices of the plane
* THREE.Face     face - Face object that denotes the sequence of vertices that make up the plane
*/
function checkPointInTriangle(point, vertices, face) {

    var linesToVertices = []

    var faceVertices = []
    faceVertices.push(vertices[face.a])
    faceVertices.push(vertices[face.b])
    faceVertices.push(vertices[face.c])

    for (var i = 0; i < faceVertices.length; i++) {
        linesToVertices.push(faceVertices[i].clone().sub(point));
    }

    var a, b, c;

    a = linesToVertices[0].angleTo(linesToVertices[1]) * (180 / Math.PI)
    b = linesToVertices[1].angleTo(linesToVertices[2]) * (180 / Math.PI)
    c = linesToVertices[2].angleTo(linesToVertices[0]) * (180 / Math.PI)

    var sum = a + b + c;

    var epsilon = 0.1;

    if (Math.abs(Math.round(sum) - 360) < epsilon) {
        return true;
    }

    return false;


}

/*
* Gets the normal given a list of vertices 
*/
function getNormal(vertices) {
    var v1 = vertices[0];
    var v2 = vertices[1];
    var v3 = vertices[2];

    var normal = (v2.clone().sub(v1).cross(v3.clone().sub(v1))).normalize();
    return normal;

}

var scene = new THREE.Scene();

var material1 = new THREE.MeshBasicMaterial({ color: 0x2194ce });
var material2 = new THREE.MeshBasicMaterial({ color: 0x4FF5ff });
var material3 = new THREE.MeshBasicMaterial({ color: 0xff0055 });

//camera
var camera = new THREE.PerspectiveCamera(75, 1, 0.1, 1000);

//light
var light = new THREE.PointLight(0xffffff, 1, 100);
light.position.set(1, 1, 4);

//plane
var planeGeometry = new THREE.PlaneGeometry(4, 20);
var plane = new THREE.Mesh(planeGeometry, material1);
plane.geometry.applyMatrix(new THREE.Matrix4().makeRotationX(-Math.PI / 2));

//sphere 1
var sphere1Geometry = new THREE.SphereGeometry(0.5, 50, 50);
var sphere1 = new THREE.Mesh(sphere1Geometry, material2);

//sphere 2
var sphere2Geometry = new THREE.SphereGeometry(0.45, 50, 50);
var sphere2 = new THREE.Mesh(sphere2Geometry, material3);


var iframe = new ImageFrame(400, 400, new THREE.Vector3(0, 0, 4));


scene.add(plane);
scene.add(sphere1);
scene.add(sphere2);
scene.add(light);

sphere1.position.set(0, 0.1, 2);
sphere2.position.set(-0.75, -0.3, 1);
plane.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(-0.5, -1, 0));
camera.position.x = 0;
camera.position.z = 6;
camera.position.y = 0;


function renderImage1(canvasId){
    
    var test = new RayTracer(camera, iframe, scene);
    test.render(canvasId);
    
}

$("#renderImage1Button").click(function(event){
    //var parentElement = $(event.toElement.parentElement);
    $(event.toElement).remove(); //remove button
    var canvasName = "renderCanvas";
    //parentElement.append("<canvas id='" + canvasName +"' style='width: 400px; height:400px'></canvas>" );
    renderImage1(canvasName);
})



//////////////////////////////////////////////////
////     HERE BE TESTING CODE FOR MY SANITY
/////////////////////////////////////////////////

/* TEST POINT IN TRIANGLE */
function testPointInTriangle() {

    var testMaterial = new THREE.MeshBasicMaterial({ color: 0x2194ce });
    var testFace = { a: 0, b: 1, c: 2 }
    var testFlatTriangleVertices = [new THREE.Vector3(0, 0, 0), new THREE.Vector3(2, 0, 0), new THREE.Vector3(0, 2, 0)]
    var testSlantTriangleVertices = [new THREE.Vector3(0, 2, -2), new THREE.Vector3(2, 0, 0), new THREE.Vector3(0, -2, 2)]

    //Checks top left corner of flat triangle and center point of the slanted triangle's straight edge
    var testPoint = new THREE.Vector3(0, 0, 0)
    
    //Check top left corner of slanted triangle
    var testPoint2 = new THREE.Vector3(0, 2, -2)
    
    //Check top edge of flat triangle
    var testPoint3 = new THREE.Vector3(1, 0, 0)
    
    //Check flat triangle center
    var testPoint4 = new THREE.Vector3(1, 0.5, 0)
    
    //Check slant triangle center
    var testPoint5 = new THREE.Vector3(0.5, 0, 0)
    
    //Check outside of bounds
    var testPoint6 = new THREE.Vector3(-1, 0, 0)
    
    var test1, test2, test3, test4, test5, test6, test7, test8;
    //top left corner flat 
    //EXPECTED TO FAIL- angle between origin vector and any other vector is NaN
    test1 = checkPointInTriangle(testPoint, testFlatTriangleVertices, testFace);
    //center long slant
    test2 = checkPointInTriangle(testPoint, testSlantTriangleVertices, testFace);
    //top left corner slant
    //EXPECTED TO FAIL - angle between origin vector and any other vector is NaN
    test3 = checkPointInTriangle(testPoint2, testSlantTriangleVertices, testFace);
    //top edge flat
    test4 = checkPointInTriangle(testPoint3, testFlatTriangleVertices, testFace);
    //flat center
    test5 = checkPointInTriangle(testPoint4, testFlatTriangleVertices, testFace);
    //slant center
    test6 = checkPointInTriangle(testPoint5, testSlantTriangleVertices, testFace);
    //out of bounds
    test7 = !checkPointInTriangle(testPoint6, testFlatTriangleVertices, testFace);
    test8 = !checkPointInTriangle(testPoint6, testSlantTriangleVertices, testFace);

    var count = 0;
    var tests = [test1,test2,test3,test4,test5,test6,test7,test8];
    
    for(var i =0; i < tests.length; i++){
        if(tests[i] == true){
            count += 1;
        } else {
            console.log("Test " + (i + 1) + " failed");
        }
    }

    console.log("Tests Succeeded: " + count + " out of 8 tests");



}

/* TEST POINT IN POLYGON */
function testPointInPolygon() {

    var testFlatPlane = {
        geometry: {
            vertices: [
                new THREE.Vector3(-1,1,0),
                new THREE.Vector3(1,1,0),
                new THREE.Vector3(-1,-1,0),
                new THREE.Vector3(1,-1,0)
            ],
            faces: [{
                a:0,
                b:1,
                c:2
            },
            {
                a:1,
                b:2,
                c:3
            }]
        }
    }
    var testSlantedPlane = {
        geometry: {
            vertices: [
                new THREE.Vector3(-1,1,-1),
                new THREE.Vector3(1,1,-1),
                new THREE.Vector3(-1,-1,1),
                new THREE.Vector3(1,-1,1)
            ],
            faces: [{
                a:0,
                b:1,
                c:2
            },
            {
                a:1,
                b:2,
                c:3
            }]
        }
    }
    
    var testPoint = new THREE.Vector3(0,0,0) //center
    var testPoint2 = new THREE.Vector3(0,1,0) // top edge flat
    var testPoint3 = new THREE.Vector3(0,1,-1) // top edge slant
    var testPoint4 = new THREE.Vector3(1,1,0) // right corner flat
    var testPoint5 = new THREE.Vector3(1,1,1) // right corner slant
    var testPoint6 = new THREE.Vector3(-2,0,0) // out of bounds
    var testPoint7 = new THREE.Vector3(0,-1,0) // bottom edge flat
    var testPoint8 = new THREE.Vector3(0,-1,1) // bottom edge slanted
    
    var test1, test2, test3, test4, test5, test6, test7, test8, test9, test10;
    
    //center
    test1 = validatePointIsInPolygonPlane(testPoint, testFlatPlane);
    test2 = validatePointIsInPolygonPlane(testPoint, testSlantedPlane);
    //top flat
    test3 = validatePointIsInPolygonPlane(testPoint2, testFlatPlane);
    //top slant
    test4 = validatePointIsInPolygonPlane(testPoint3, testSlantedPlane);
    //corner flat
    //EXPECTED TO FAIL- angle between origin vector and any other vector is NaN
    test5 = validatePointIsInPolygonPlane(testPoint4, testFlatPlane);
    //corner slant
    //EXPECTED TO FAIL- angle between origin vector and any other vector is NaN
    test6 = validatePointIsInPolygonPlane(testPoint5, testSlantedPlane);
    //out of bounds
    test7 = !validatePointIsInPolygonPlane(testPoint6, testFlatPlane);
    test8 = !validatePointIsInPolygonPlane(testPoint6, testSlantedPlane);
    
    //bottom edges
    test9 = validatePointIsInPolygonPlane(testPoint7, testFlatPlane);
    test10 = validatePointIsInPolygonPlane(testPoint8, testSlantedPlane);

    var count = 0;
    var tests = [test1,test2,test3,test4,test5,test6,test7,test8, test9, test10];
    
    for(var i =0; i < tests.length; i++){
        if(tests[i] == true){
            count += 1;
        } else {
            console.log("Test " + (i + 1) + " failed");
        }
    }
    
    console.log("Tests Succeeded: " + count + " out of " + tests.length + " tests");

}

/* TEST SPHERE INTERSECTION */
function testSphereIntersection() {

    var testMaterial = new THREE.MeshBasicMaterial({ color: 0x2194ce });
    var testSphereGeometry = new THREE.SphereGeometry(0.5, 50, 50);
    
    //1 unit diameter sphere @ origin
    var testSphere = new THREE.Mesh(testSphereGeometry, testMaterial);
    
    //ray looking at origin, looking 5 Z-units away. 4 units away from closest point on sphere
    var testRay = new Ray(
        new THREE.Vector3(0, 0, 5),
        new THREE.Vector3(0, 0, -5).normalize()
        )
    var expectedIntersectionForTestRay = new THREE.Vector3(0, 0, 0.5);
    
    //ray looking at right most part of the sphere
    var testRay2 = new Ray(
        new THREE.Vector3(0, 0, 5),
        new THREE.Vector3(0.5, 0, -5).normalize()
        )
    var expectedIntersectionForTestRay2 = new THREE.Vector3(0.5, 0, 0);
    //ray looking at right most part of the sphere
    var testRay3 = new Ray(
        new THREE.Vector3(0, 0, 5),
        new THREE.Vector3(1, 0, -5).normalize()
        )

    var test1 = getSphereIntersection(testRay, testSphere);
    var test2 = getSphereIntersection(testRay2, testSphere);
    var test3 = getSphereIntersection(testRay3, testSphere);

    if (test1 != null && test2 != null) {

        console.log("Test 1 Intersection: (" + test1.x + "," + test1.y + "," + test1.z + ") \n Expected: "
            + " (" + expectedIntersectionForTestRay.x + "," + expectedIntersectionForTestRay.y + "," + expectedIntersectionForTestRay.z + ")")

        console.log("Test 2 Intersection: (" + test2.x + "," + test2.y + "," + test2.z + ") \n Expected: "
            + " (" + expectedIntersectionForTestRay2.x + "," + expectedIntersectionForTestRay2.y + ", ~" + expectedIntersectionForTestRay2.z + ")")

    }
    
    if(test3 != null){
        console.log("Error: Test 3 found an intersection point");
    }

}

/* TEST PLANE INTERSECTION */
function testPlaneIntersection() {
    
    var testFlatPlane = {
        geometry: {
            vertices: [
                new THREE.Vector3(-1,1,0),
                new THREE.Vector3(-1,-1,0),
                new THREE.Vector3(1,1,0),
                new THREE.Vector3(1,-1,0)
            ],
            faces: [{
                a:0,
                b:1,
                c:2
            },
            {
                a:1,
                b:2,
                c:3
            }]
        }
    }
    var testSlantedPlane = {
        geometry: {
            vertices: [
                new THREE.Vector3(-1,1,-1),
                new THREE.Vector3(-1,-1,1),
                new THREE.Vector3(1,1,-1),
                new THREE.Vector3(1,-1,1)
            ],
            faces: [{
                a:0,
                b:1,
                c:2
            },
            {
                a:1,
                b:2,
                c:3
            }]
        }
    }
    
    //center
    var testRay = new Ray(
        new THREE.Vector3(0,0,5),
        new THREE.Vector3(0,0,-5).normalize()
    )
    
    //flat top
    var testRay2 = new Ray(
        new THREE.Vector3(0,0,5),
        new THREE.Vector3(0,1,-5).normalize()
    )
    
    //slant top
    var testRay3 = new Ray(
        new THREE.Vector3(0,0,5),
        new THREE.Vector3(0,1,-6).normalize()
    )

    //out of bounds
    var testRay4 = new Ray(
        new THREE.Vector3(0,0,5),
        new THREE.Vector3(0,5, 0).normalize()
    )
    
    var test1,test2,test3,test4,test5,test6,test7,test8;
    var test9,test10,test11,test12;
    
    //center
    test1 = getPlaneIntersection(testRay, testFlatPlane);
    test2 = getPlaneIntersection(testRay, testSlantedPlane);
    
    test3 = getPlaneIntersectionNotSchool(testRay, testFlatPlane);
    test4 = getPlaneIntersectionNotSchool(testRay, testSlantedPlane);
    
    
    //flat top
    test5 = getPlaneIntersection(testRay2, testFlatPlane);
    test6 = getPlaneIntersection(testRay3, testSlantedPlane);
    
    //slanted top
    test7 = getPlaneIntersectionNotSchool(testRay2, testFlatPlane);
    test8 = getPlaneIntersectionNotSchool(testRay3, testSlantedPlane);
    
    //out of bounds
    
    test9 = getPlaneIntersection(testRay4, testFlatPlane);
    test10 = getPlaneIntersection(testRay4, testSlantedPlane);
    
    test11 = getPlaneIntersectionNotSchool(testRay4, testFlatPlane);
    test12 = getPlaneIntersectionNotSchool(testRay4, testSlantedPlane);

    
    var count = 0;
    var tests = [test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,test11,test12];
    
    for(var i =0; i < tests.length; i++){
          
        if(tests[i] != null){    
            if(i % 2 == 0 && validatePlaneEquation(getNormal(testFlatPlane.geometry.vertices.slice(0,4)), testFlatPlane.geometry.vertices[0], tests[i])){
                count += 1;
                
            } else if(validatePlaneEquation(getNormal(testSlantedPlane.geometry.vertices.slice(0,4)), testSlantedPlane.geometry.vertices[0], tests[i])){
                count += 1;
            } else {
                console.log("Test " + (i + 1) + " failed");
            }
        } else {
            console.log("Test " + (i + 1) + " failed");
        }
    }
    
    console.log("Tests Succeeded: " + count + " out of " + tests.length + " tests");
    
}

function validatePlaneEquation(normal, point, intersectPoint){
    
    var d = normal.clone().dot(intersectPoint.clone().sub(point));
    var epsilon = 0.1;
    if(Math.abs(d) < epsilon){
        return true;
    }
    
    return false;
    
}