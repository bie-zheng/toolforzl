from collections import Counter
from ovito.io import import_file
from ovito.pipeline import FileSource
from ovito.modifiers import (
    ClusterAnalysisModifier,
    ExpressionSelectionModifier,
    WignerSeitzAnalysisModifier,
)

def count_interstitial_atoms(reference_file, target_file, **kwargs):
    """
    计算间隙原子数量

    参数:
        reference_file (str): 参考构型文件路径
        target_file (str): 目标构型文件路径
        **kwargs: 传递给 WignerSeitzAnalysisModifier 的额外关键字参数
                  (默认 output_displaced=False，可通过此参数覆盖)

    返回:
        tuple: (间隙原子数, 空位数)
    """
    # 同时加载参考文件和目标文件作为两帧（frame 0 = 参考, frame 1 = 目标）
    ref_source = FileSource()
    ref_source.load(reference_file)
    pipeline = import_file(target_file)
    

    # 创建并应用Wigner-Seitz缺陷分析修改器
    ws_params = {'output_displaced': False, **kwargs}
    ws_modifier = WignerSeitzAnalysisModifier(**ws_params)
    ws_modifier.reference = ref_source
    pipeline.modifiers.append(ws_modifier)

    # 计算第1帧（目标文件），自动以第0帧为参考
    data = pipeline.compute(1)
    interstitial_count = data.attributes['WignerSeitz.interstitial_count']


    return interstitial_count

def count_interstitial_atoms(reference_file, target_file, **kwargs):
    """
    计算间隙原子数量

    参数:
        reference_file (str): 参考构型文件路径
        target_file (str): 目标构型文件路径
        **kwargs: 传递给 WignerSeitzAnalysisModifier 的额外关键字参数
                  (默认 output_displaced=False，可通过此参数覆盖)

    返回:
        tuple: (间隙原子数, 空位数)
    """
    # 同时加载参考文件和目标文件作为两帧（frame 0 = 参考, frame 1 = 目标）
    ref_source = FileSource()
    ref_source.load(reference_file)
    pipeline = import_file(target_file)
    

    # 创建并应用Wigner-Seitz缺陷分析修改器
    ws_params = {'output_displaced': False, **kwargs}
    ws_modifier = WignerSeitzAnalysisModifier(**ws_params)
    ws_modifier.reference = ref_source
    pipeline.modifiers.append(ws_modifier)

    # 计算第1帧（目标文件），自动以第0帧为参考
    data = pipeline.compute(1)    
    vacancy_count = data.attributes['WignerSeitz.vacancy_count']

    return vacancy_count


def count_interstitial_clusters(reference_file, target_file, cutoff=3.0):
    """
    计算目标文件相对于参考帧的间隙原子团簇大小分布。
    WS分析使用atom，并且包括了在间隙原子位点上的原位原子。（这对于计算辐照损伤的后果是合理的）

    参数:
        reference_file (str): 参考构型文件路径
        target_file (str): 目标构型文件路径
        cutoff (float): 团簇截止半径（单位：Å），默认 3.0

    返回:
        dict: {团簇大小: 该大小的团簇数量}，如 {1: 5, 2: 3}
    """
    ref_source = FileSource()
    ref_source.load(reference_file)

    pipeline = import_file(target_file)
    # 使用 Wigner-Seitz 缺陷分析识别间隙原子
    ws_modifier = WignerSeitzAnalysisModifier(output_displaced=True)
    #output_display=1用于设置计算的是atom类型
    ws_modifier.reference = ref_source
    pipeline.modifiers.append(ws_modifier)

    # 选中所有间隙原子（Occupancy >= 2）
    pipeline.modifiers.append(
        ExpressionSelectionModifier(expression="Occupancy >= 2")
    )

    # 仅对选中的间隙原子进行团簇分析
    cluster_modifier = ClusterAnalysisModifier(
        cutoff=cutoff,
        only_selected=True,
    )
    pipeline.modifiers.append(cluster_modifier)

    data = pipeline.compute(1)
    cluster_ids = data.particles["Cluster"]

    # 统计每个团簇包含的原子数
    cluster_sizes = Counter(cluster_ids)

    # 移除未聚类的噪声点（Cluster ID 为 0）
    cluster_sizes.pop(0, None)

    # 统计每种大小的团簇数量
    size_distribution = Counter(cluster_sizes.values())
    return dict(sorted(size_distribution.items()))


def count_vacancy_clusters(reference_file, target_file, cutoff=3.0):
    """
    计算目标文件相对于参考帧的空位团簇大小分布。
    WS分析使用 lattice site 位点，选择 Occupancy == 0 的空位位置。

    参数:
        reference_file (str): 参考构型文件路径
        target_file (str): 目标构型文件路径
        cutoff (float): 团簇截止半径（单位：Å），默认 3.0

    返回:
        dict: {团簇大小: 该大小的团簇数量}，如 {1: 5, 2: 3}
    """
    ref_source = FileSource()
    ref_source.load(reference_file)

    pipeline = import_file(target_file)
    # 使用 Wigner-Seitz 缺陷分析识别空位，output_displaced=False 输出晶格位点
    ws_modifier = WignerSeitzAnalysisModifier(output_displaced=False)
    ws_modifier.reference = ref_source
    pipeline.modifiers.append(ws_modifier)

    # 选中所有空位（Occupancy == 0）
    pipeline.modifiers.append(
        ExpressionSelectionModifier(expression="Occupancy == 0")
    )

    # 仅对选中的空位进行团簇分析
    cluster_modifier = ClusterAnalysisModifier(
        cutoff=cutoff,
        only_selected=True,
    )
    pipeline.modifiers.append(cluster_modifier)

    data = pipeline.compute(1)
    cluster_ids = data.particles["Cluster"]

    # 统计每个团簇包含的位数
    cluster_sizes = Counter(cluster_ids)

    # 移除未聚类的噪声点（Cluster ID 为 0）
    cluster_sizes.pop(0, None)

    # 统计每种大小的团簇数量
    size_distribution = Counter(cluster_sizes.values())
    return dict(sorted(size_distribution.items()))


# 使用示例
if __name__ == "__main__":
    reference = fr"F:\科研\202607\FeCrAl_bcc_sphere.data"  # 替换为实际参考文件路径
    target = fr"F:\科研\202607\600.38597.10000.2500.dump"  # 替换为实际目标文件路径
    
    interstitials, vacancies = count_interstitial_atoms(reference, target)
    print(f"间隙原子数: {interstitials}, 空位数: {vacancies}")

    reference = r"F:\科研\202607\600.38597.10000.0.dump"
    target = r"F:\科研\202607\600.38597.10000.2500.dump"

    size_distribution = count_interstitial_clusters(reference, target, cutoff=3.0)
    print(f"间隙原子团簇大小分布: {size_distribution}")

    reference = r"F:\科研\202607\600.38597.10000.0.dump"
    target = r"F:\科研\202607\600.38597.10000.2500.dump"

    vacancy_size_distribution = count_vacancy_clusters(reference, target, cutoff=3.0)
    print(f"空位团簇大小分布: {vacancy_size_distribution}")