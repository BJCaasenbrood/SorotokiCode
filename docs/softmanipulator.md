# Soft manipulator

<div align="center">
<video width="400" height="310" autoplay controls loop>
<source src="assets/fem_crawler.mp4" type="video/mp4">
</video>
</div>

<script type="text/javascript" src="https://static.sketchfab.com/api/sketchfab-viewer-1.12.1.js"></script>

<div align="center">
<body>
    <!-- Insert an empty iframe with attributes -->
    <iframe width="850" height="750"src="" id="api-frame" allow="autoplay; fullscreen; xr-spatial-tracking" xr-spatial-tracking execution-while-out-of-viewport execution-while-not-rendered web-share allowfullscreen mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

    <!-- Initialize the viewer -->
    <script type="text/javascript">
    var iframe = document.getElementById( 'api-frame' );
    var uid = '3fd62289149f4145856b9c664c033917';

    // By default, the latest version of the viewer API will be used.
    var client = new Sketchfab( iframe );

    // Alternatively, you can request a specific version.
    // var client = new Sketchfab( '1.12.1', iframe );

    client.init( uid, {
        success: function onSuccess( api ){
            api.start();
            api.addEventListener( 'viewerready', function() {

                // API is ready to use
                // Insert your code here
                console.log( 'Viewer is ready' );

            } );
        },
        error: function onError() {
            console.log( 'Viewer error' );
        }
    } );
    </script>

</body>
</div>
