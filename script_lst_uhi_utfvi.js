// 
// Análisis de Temperatura Superficial (LST), Isla de Calor Urbano (UHI)
// e Índice de Variación del Campo Térmico Urbano (UTFVI)
//
// Región Metropolitana de Santiago, Chile.
// 
// Satélite: Landsat 8 Collection 2
// 
// Autor: Lucas Blachet Del Solar
// Fecha: Noviembre 2025
// Email: lucas.blachet@uc.cl
//

// INPUTS
var ae = ee.FeatureCollection('projects/dynamic-world-test-1/assets/CL_RM');
var fecha_inicio = '2024-12-01';
var fecha_fin = '2025-03-31';

Map.centerObject(ae, 9);

// FUNCIONES

// Factores de escala USGS Landsat 8 Collection 2
function FactoresEscala(image) {
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
    return image.addBands(opticalBands, null, true)
        .addBands(thermalBands, null, true);
}

// Máscara de nubes usando QA_PIXEL
function MascaraNubes(col) {
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    var qa = col.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    return col.updateMask(mask);
}

// PROCESAMIENTO
var image = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(fecha_inicio, fecha_fin)
    .filterBounds(ae)
    .map(FactoresEscala)
    .map(MascaraNubes)
    .median()
    .clip(ae);

// Visualización RGB
var visualization = {
    bands: ['SR_B4', 'SR_B3', 'SR_B2'],
    min: 0.0,
    max: 0.3,
};
Map.addLayer(image, visualization, 'RGB', false);

// NDVI
var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
Map.addLayer(ndvi, { min: -1, max: 1, palette: ['blue', 'white', 'green'] }, 'NDVI', false);

// Estadísticas NDVI
var ndvi_minimo = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: ae,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

var ndvi_maximo = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: ae,
    scale: 30,
    maxPixels: 1e9
}).values().get(0));

// Fracción de vegetación
var fv = (ndvi.subtract(ndvi_minimo).divide(ndvi_maximo.subtract(ndvi_minimo)))
    .pow(ee.Number(2))
    .rename('FV');

// Emisividad
var em = fv.multiply(ee.Number(0.004)).add(ee.Number(0.986)).rename('EM');

// Temperatura de brillo
var thermal = image.select('ST_B10').rename('thermal');

// LST (Land Surface Temperature)
var lst = thermal.expression(
    '(tb / (1 + ((11.5 * (tb / 14380)) * log(em)))) - 273.15',
    {
        'tb': thermal.select('thermal'),
        'em': em
    }
).rename('LST');

// Paleta térmica
var lst_vis = {
    min: 10,
    max: 45,
    palette: [
        '313695', '4575b4', '74add1', 'abd9e9', 'e0f3f8',
        'ffffbf', 'fee090', 'fdae61', 'f46d43', 'd73027', 'a50026'
    ]
};

Map.addLayer(lst, lst_vis, 'LST (°C)');

// ESTADÍSTICAS LST

var stats_lst = lst.reduceRegion({
    reducer: ee.Reducer.mean()
        .combine(ee.Reducer.stdDev(), '', true)
        .combine(ee.Reducer.min(), '', true)
        .combine(ee.Reducer.max(), '', true)
        .combine(ee.Reducer.percentile([25, 50, 75]), '', true),
    geometry: ae,
    scale: 30,
    maxPixels: 1e9
});

var lst_mean = ee.Number(stats_lst.get('LST_mean'));
var lst_std = ee.Number(stats_lst.get('LST_stdDev'));
var lst_min = ee.Number(stats_lst.get('LST_min'));
var lst_max = ee.Number(stats_lst.get('LST_max'));
var lst_p25 = ee.Number(stats_lst.get('LST_p25'));
var lst_p50 = ee.Number(stats_lst.get('LST_p50'));
var lst_p75 = ee.Number(stats_lst.get('LST_p75'));

print('ESTADÍSTICAS LST - SANTIAGO');
print('Media:', lst_mean, '°C');
print('Desv. Estándar:', lst_std, '°C');
print('Mínima:', lst_min, '°C');
print('Máxima:', lst_max, '°C');
print('Percentil 25:', lst_p25, '°C');
print('Mediana (P50):', lst_p50, '°C');
print('Percentil 75:', lst_p75, '°C');

// UHI (Urban Heat Island)
var uhi = lst.subtract(lst_mean).divide(lst_std).rename('UHI');

var uhi_vis = {
    min: -3,
    max: 3,
    palette: ['313695', '4575b4', 'abd9e9', 'ffffbf', 'fdae61', 'f46d43', 'd73027']
};

Map.addLayer(uhi, uhi_vis, 'UHI');

// UTFVI (Urban Thermal Field Variance Index)
var utfvi = lst.subtract(lst_mean)
    .divide(lst_std)
    .rename('UTFVI');

var utfvi_vis = {
    min: -3,
    max: 3,
    palette: ['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
};

Map.addLayer(utfvi, utfvi_vis, 'UTFVI');

// Exportar LST
Export.image.toDrive({
    image: lst,
    description: 'LST_RM',
    folder: 'GEE_Exports',
    region: ae,
    scale: 30,
    maxPixels: 1e13,
    crs: 'EPSG:4326'
});

// Exportar UHI
Export.image.toDrive({
    image: uhi,
    description: 'UHI_RM',
    folder: 'GEE_Exports',
    region: ae,
    scale: 30,
    maxPixels: 1e13,
    crs: 'EPSG:4326',
});

// Exportar UTFVI
Export.image.toDrive({
    image: utfvi,
    description: 'UTFVI_RM',
    folder: 'GEE_Exports',
    region: ae,
    scale: 30,
    maxPixels: 1e13,
    crs: 'EPSG:4326'
});

// GRÁFICOS Y VISUALIZACIONES

// Histograma de LST
var histogram_lst = ui.Chart.image.histogram({
    image: lst,
    region: ae,
    scale: 30,
    maxPixels: 1e13
}).setOptions({
    title: 'Distribución de Temperatura Superficial (LST)',
    hAxis: { title: 'Temperatura (°C)' },
    vAxis: { title: 'Frecuencia (píxeles)' },
    colors: ['#d73027'],
    legend: { position: 'none' }
});
print(histogram_lst);

// Histograma de UHI
var histogram_uhi = ui.Chart.image.histogram({
    image: uhi,
    region: ae,
    scale: 30,
    maxPixels: 1e13
}).setOptions({
    title: 'Distribución de Isla de Calor Urbana (UHI)',
    hAxis: { title: 'UHI (desviaciones estándar)' },
    vAxis: { title: 'Frecuencia (píxeles)' },
    colors: ['#4575b4'],
    legend: { position: 'none' }
});
print(histogram_uhi);

// Relación NDVI vs LST (scatter plot con muestreo)
var sample = ndvi.addBands(lst).sample({
    region: ae,
    scale: 120,
    numPixels: 5000,
    seed: 42,
    geometries: true
});

var scatter_ndvi_lst = ui.Chart.feature.byFeature(sample, 'NDVI', 'LST')
    .setChartType('ScatterChart')
    .setOptions({
        title: 'Relación NDVI vs LST',
        hAxis: { title: 'NDVI' },
        vAxis: { title: 'LST (°C)' },
        pointSize: 2,
        colors: ['#238b45'],
        trendlines: { 0: { color: 'red', lineWidth: 2, opacity: 0.5 } }
    });
print(scatter_ndvi_lst);
